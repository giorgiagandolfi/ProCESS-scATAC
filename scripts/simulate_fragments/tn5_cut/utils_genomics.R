library(data.table)
binRegion = function(start, end, binSize=NULL, binCount=NULL, indicator=NULL) {
  if (is.null(binSize) & is.null(binCount)) {
    stop("You must provide either binSize or binCount")
  }
  if (is.null(binSize)) {
    binSize = round(sum(end-start)/binCount)
  }
  binCountByChrom = round((end-start)/binSize)
  binCountByChrom[binCountByChrom==0]=1
  binSizeByChrom = (end-start)/(binCountByChrom)
  breaks = round(unlist(lapply(binCountByChrom, 
                               function(x) seq(from=0, to=x))) * 
                   rep(binSizeByChrom, (binCountByChrom+1)))
  endpoints = cumsum(binCountByChrom + 1) 
  startpoints = c(1, endpoints[-length(endpoints)]+1)
  
  dataTable = data.table(start=breaks[-endpoints]+1, 
                         end=breaks[-startpoints],
                         id=rep((seq_along(start)), binCountByChrom),
                         binID=unlist(lapply(binCountByChrom, 
                                             function(x) seq(from=1, to=x))),
                         ubinID=seq_along(breaks[-startpoints]),
                         key="id")
  
  if (!is.null(indicator)){
    idCol = rep(indicator, binCountByChrom)
    dataTable = data.table(idCol, dataTable)
  }
  return(dataTable)
}