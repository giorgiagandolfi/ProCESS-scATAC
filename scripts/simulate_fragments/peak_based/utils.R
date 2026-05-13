library(extraDistr)

rztpois <- function(n, lambda) {
  x <- rpois(n, lambda)
  
  # resample zeros
  while(any(x == 0)) {
    idx <- which(x == 0)
    x[idx] <- rpois(length(idx), lambda)
  }
  
  x
}

#### Function to sample fragment counts for each peak
sample_fragment_count <- function(lambda = 0.1,
                                  max_count = 20) {
  
  repeat {
    
    x <- rztpois(1, lambda)
    
    if (x <= max_count) {
      return(x)
    }
  }
}


#### Functionto sample fragment sizes
sample_fragment_size <- function() {
  
  comp <- sample(
    x = 1:4,
    size = 1,
    prob = c(0.5, 0.35, 0.10, 0.05)
  )
  
  if (comp == 1) {
    size <- rtnorm(1, mean = 50, sd = 40,a = 1,b=150)
  } else if (comp == 2) {
    size <- rnorm(1, mean = 200, sd = 50)
  } else if (comp == 3) {
    size <- rnorm(1, mean = 400, sd = 30)
  } else {
    size <- rnorm(1, mean = 600, sd = 40)
  }
  
  round(max(size, 1))
}

##### Place fragments 

place_fragments <- function(start, fragments_sizes) {
  n <- length(fragments_sizes)
  from <- numeric(n)
  to <- numeric(n)
  
  current <- start
  
  for (i in seq_along(fragments_sizes)) {
    from[i] <- current
    to[i] <- current + fragments_sizes[i] - 1
    current <- to[i] + 1
  }
  fragm_ids <- paste0("fragm_",seq_along(1:length(fragments_sizes)))
  frg_pos_df <-data.frame(fragment = fragm_ids,
                          from = from,
                          to = to)
  return(frg_pos_df)
}

place_fragments_in_peak <- function(fragments_sizes, peak_id,peak_from,peak_to){
  # peak <- crc_peaks[1,]
  # peak_from <- peak$from
  # peak_to <- peak$to
  peak_center <- peak_from+((peak_to -peak_from)/2)
  
  fragments_center <- rnorm(n=1,mean = as.numeric(peak_center),sd = 10)
  sum_fragments <- sum(fragments_sizes)
  superfrag_from <- fragments_center-(sum_fragments/2)
  superfrag_to <- fragments_center+(sum_fragments/2)
  
  fragm_ids <- paste0("fragm_",seq_along(1:length(fragments_sizes)))
  frg_pos_df <- place_fragments(start = superfrag_from,fragments_sizes = fragments_sizes)
  frg_pos_df <- frg_pos_df %>% 
    mutate(peak=peak_id)
  
}





###### Conditional sampling
sample_fragments_for_peak <- function(
    #peak_size = 500,
    peak_id, ## in format chr1:242947:562757
    peak_from,
    peak_to,
    flank = 50,
    lambda = 0.1,
    max_attempts = 100,
    fragment_len_dist=NULL,
    tot_cn=2
) {
  peak_size <- peak_to-peak_from
  
  max_total <- peak_size + 2 * flank
  
  for (attempt in 1:max_attempts) {
    
    n_frag <- sample_fragment_count(lambda)
    
    
    success <- TRUE
    allele_size_list=list()
    for (allele in 1:tot_cn){
      sizes <- numeric(0)
      for (i in 1:n_frag) {
        
        accepted <- FALSE
        
        for (j in 1:max_attempts) {
          
          # 
          if (!is.null(fragment_len_dist)){
            s <- sample(x = fragment_len_dist,size = 1)
          } else{
            s <- sample_fragment_size()
          }
          
          # reserve at least 1bp for remaining fragments
          min_needed <- (n_frag - i)
          
          if ((sum(sizes) + s + min_needed) <= max_total) {
            sizes <- c(sizes, s)
            accepted <- TRUE
            break
          }
        }
        
        if (!accepted) {
          success <- FALSE
          break
        }
      }
      
      if (success) {
        
        peaks_frags_df <- place_fragments_in_peak(fragments_sizes = sizes,
                                                  peak_id = peak_id,peak_from = peak_from,
                                                  peak_to = peak_to)
        allele_size_list[[allele]]=peaks_frags_df %>% mutate(allele=paste0('allele_',allele))
        # return(peaks_frags_df)
        
      }
    }

  }
  results_allele=do.call('rbind',allele_size_list)
  return(results_allele)
  # stop("Could not generate valid configuration")
}


