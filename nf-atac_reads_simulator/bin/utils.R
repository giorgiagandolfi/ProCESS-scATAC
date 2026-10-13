# library(extraDistr)
library(data.table)
library(Matrix)
library(patchwork)
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
                                  max_count = 20,accessibility_score=NULL) {
  
  repeat {
    
    x <- rztpois(1, (lambda*accessibility_score))
    
    if (x <= max_count) {
      return(x)
    }
  }
}


#### Functionto sample fragment sizes
sample_fragment_size <- function(n) {
  
  comp <- sample(
    x = 1:4,
    size = n, # change to n fragments
    prob = c(0.5, 0.35, 0.10, 0.05)
  )
  
  if (comp == 1) {
    size <- rnorm(1, mean = 50, sd = 40,a = 1,b=150)
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
### old version
place_fragments_old <- function(start, fragments_sizes,gaps_sizes) {
  n <- length(fragments_sizes)
  from <- numeric(n)
  to <- numeric(n)
  
  current <- start
  if (n>1){
    gaps_sizes=gaps_sizes[-1]
  }
  for (i in seq_along(fragments_sizes)) {
    if (i == 1){
      from[i] <- current
      to[i] <- current + fragments_sizes[i] - 1
      if (from[i]>=to[i]){
        message("Error")
      }
      current <- to[i] + 1
    } else{
      for (j in seq_along(gaps_sizes)){
        from[i] <- current + gaps_sizes[j]
        to[i] <- current + fragments_sizes[i] - 1
        current <- to[i] + 1
      }
    }
  }
  fragm_ids <- paste0("fragm_",seq_along(1:length(fragments_sizes)))
  frg_pos_df <-data.frame(fragment = fragm_ids,
                          from = from,
                          to = to)
  return(frg_pos_df)
}

place_fragments <- function(start, fragments_sizes, gaps_sizes) {
  
  n <- length(fragments_sizes)
  
  if (length(gaps_sizes) != n) {
    stop("gaps_sizes must have same length as fragments_sizes")
  }
  
  from <- numeric(n)
  to <- numeric(n)
  
  from[1] <- start
  to[1] <- start + fragments_sizes[1] - 1
  
  if (n > 1) {
    for (i in 2:n) {
      from[i] <- to[i - 1] + gaps_sizes[i] + 1
      to[i] <- from[i] + fragments_sizes[i] - 1
    }
  }
  
  data.frame(
    fragment = paste0("fragm_", seq_len(n)),
    from = from,
    to = to
  )
}

place_fragments_in_peak <- function(fragments_sizes, gaps_sizes,peak_id,peak_from,peak_to,sd_peak_center=30){
  # peak <- crc_peaks[1,]
  # peak_from <- peak$from
  # peak_to <- peak$to
  peak_center <- peak_from+((peak_to -peak_from)/2)
  
  fragments_center <- rnorm(n=1,mean = as.numeric(peak_center),sd = sd_peak_center)
  sum_fragments <- sum(fragments_sizes,gaps_sizes)
  superfrag_from <- fragments_center-(sum_fragments/2)
  superfrag_to <- fragments_center+(sum_fragments/2)
  
  fragm_ids <- paste0("fragm_",seq_along(1:length(fragments_sizes)))
  frg_pos_df <- place_fragments(start = superfrag_from,fragments_sizes = fragments_sizes,gaps_sizes=gaps_sizes)
  frg_pos_df <- frg_pos_df %>% 
    mutate(peak=peak_id)
  
}

sample_gap_size <- function(mu = 50, sigma = 10) {
  g <- rnorm(1, mean = mu, sd = sigma)
  g <- round(g) 
  g <- max(0, g)
  return(g)
}



###### Conditional sampling
sample_fragments_for_peak <- function(
    #peak_size = 500,
    peak_id, ## in format chr1:242947:562757
    peak_from,
    peak_to,
    flank = 100,
    lambda = 0.2, #0.1,
    max_attempts = 100,
    fragment_len_dist=NULL,
    tot_cn=2
    # peak_status="open"
) {
    peak_size <- peak_to-peak_from
    
    max_total <- peak_size + 2 * flank
    
    for (attempt in 1:max_attempts) {
      
      n_frag <- sample_fragment_count(lambda)
      # n_frag - 2
      
      success <- TRUE
      allele_size_list=list()
      for (allele in 1:tot_cn){
        sizes <- numeric(0)
        gaps <- numeric(0)
        for (i in 1:n_frag) {
          
          accepted <- FALSE
          
          for (j in 1:max_attempts) {
            
            # 
            if (!is.null(fragment_len_dist)){
              s <- sample(x = fragment_len_dist,size = 1)
            } else{
              s <- sample_fragment_size()
            }
            
            # sample gap BEFORE fragment (except first fragment)
            g <- 0
            if (i > 1) {
              g <- sample_gap_size()
            }
            
            remaining_frags <- n_frag - i
            remaining_gaps  <- max(0, n_frag - i - 1)
            
            min_needed <- remaining_frags + remaining_gaps  # or more realistic below
            
            if ((sum(sizes) + sum(gaps) + s + g + min_needed) <= max_total) {
              
              sizes <- c(sizes, s)
              if (i > 1) {
                gaps <- c(gaps, g)
              }
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
                                                    gaps_sizes = gaps,
                                                    peak_id = peak_id,peak_from = peak_from,
                                                    peak_to = peak_to)
          allele_size_list[[allele]]=peaks_frags_df %>% mutate(allele=paste0('allele_',allele))
        }
      }
      
  }
  results_allele=do.call('rbind',allele_size_list)
  return(results_allele)
  # stop("Could not generate valid configuration")
}

sample_fragments_for_peak_vec <- function(
    peak_id,
    peak_chr,
    peak_from,
    peak_to,
    flank = 100,
    lambda = 0.2,
    fragment_len_dist = NULL,
    tot_cn,
    max_fragments = 50,
    cell_id
) {
  
  n <- length(peak_id)
  
  stopifnot(
    length(peak_from) == n,
    length(peak_to) == n,
    length(tot_cn) == n
  )
  
  peak_size <- peak_to - peak_from
  max_total <- peak_size + 2 * flank
  
  # vectorized fragment count sampling
  n_frag_vec <- rztpois(n = n,lambda = lambda)
  
  n_frag_vec <- pmin(n_frag_vec, max_fragments)
  
  results <- vector("list", n)
  
  for (k in seq_len(n)) {
    n_frag <- n_frag_vec[k]
    
    if (n_frag == 0 || tot_cn[k] == 0) {
      next
    }
    
    allele_results <- vector("list", tot_cn[k])
    
    for (allele in seq_len(tot_cn[k])) {
      
      # sample all fragment sizes simultaneously
      if (!is.null(fragment_len_dist)) {
        
        sizes <- sample(
          fragment_len_dist,
          size = n_frag,
          replace = TRUE
        )
        
      } else {
        
        sizes <- sample_fragment_size(n_frag)
        
      }
      
      # sample all gaps simultaneously
      gaps <- numeric(n_frag)
      
      if (n_frag > 1) {
        gaps[2:n_frag] <- sample_gap_size(n_frag - 1)
      }
      
      # cumulative occupied space
      occupied <- cumsum(sizes + gaps)
      
      keep <- occupied <= max_total[k]
      
      if (!any(keep)) {
        next
      }
      
      sizes <- sizes[keep]
      gaps  <- gaps[keep]
      print(peak_id[k])
      print(gaps)
      print(sizes)
      # vectorized fragment placement
      peaks_frags_df <- place_fragments_in_peak(fragments_sizes = sizes,
                                                gaps_sizes = gaps,
                                                peak_id = peak_id[k],peak_from = peak_from[k],
                                                peak_to = peak_to[k],sd_peak_center=100)
      
      allele_results[[allele]] <- peaks_frags_df %>% mutate(fragment_allele=paste0('allele_',allele)) %>% 
        mutate(peak_id = peak_id[k],
               peak_chr=peak_chr,
               peak_from = peak_from[k],
               peak_to = peak_to[k],
               fragment_size = sizes,
               cell_id=cell_id,
               fragment_chr=peak_chr
               ) %>% 
        dplyr::rename(fragment_start=from,
                      fragment_end=to)

    }
    
    results[[k]] <- do.call("rbind",allele_results)
  }
  
  do.call("rbind",results)
}


# get_background_regions <- function(peak_fragments_df,genome,gaps_file,centromeres_file,filter_small_than=150){
#   #### peak_df is a dataframe with all the peaks for that cell
#   #### like this one
#   # peak_id peak_chr peak_start  peak_end
#   # 1        chr1-1006219-1006719     chr1    1006219   1006719
#   # 2        chr1-1068981-1069481     chr1    1068981   1069481
#   # 3      chr1-10694533-10695033     chr1   10694533  10695033
#   # 4    chr1-108559703-108560203     chr1  108559703 108560203
#   # 5    chr1-109763665-109764165     chr1  109763665 109764165
#   # 6    chr1-111204245-111204745     chr1  111204245 111204745
#   peak_region_df <- peak_fragments_df %>% 
#     group_by(peak) %>% 
#     summarise(start_peak_region=min(fragment_start),
#               end_peak_region=max(fragment_end),
#               chr_peak_region=paste0('chr',unique(fragment_chr)))
#   gr_peak_union = GRanges(
#     seqnames = peak_region_df$chr_peak_region,
#     ranges = IRanges(start = peak_region_df$start_peak_region, end = peak_region_df$end_peak_region),
#   )
#   if (genome=='hg38'){
#     gr_genome <- GRanges(
#       seqnames = names(seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)),
#       ranges = IRanges(
#         start = 1,
#         end = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
#       ))
#   } else if (genome=='hg19'){
#     gr_genome <- GRanges(
#       seqnames = names(seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)),
#       ranges = IRanges(
#         start = 1,
#         end = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
#       ))
#   }
#   gr_outside_peaks <- setdiff(gr_genome, gr_peak_union)
#   
#   blacklist_gr <- BiocIO::import(
#     "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
#   )
#   
#   blacklist_gr <- keepStandardChromosomes(
#     blacklist_gr,
#     pruning.mode = "coarse"
#   )
#   
#   cytoband <- read.table(
#     gzfile(centromeres_file), ## must be a gz
#     header = FALSE,
#     sep = '\t',
#     col.names = c("chr", "start", "end", "name", "stain")
#   )
#   centromeres <- GRanges(
#     seqnames = cytoband$chr[cytoband$stain == "acen"],
#     ranges = IRanges(
#       start = cytoband$start[cytoband$stain == "acen"] + 1, # BED is 0-based
#       end = cytoband$end[cytoband$stain == "acen"]
#     )
#   )
#   
#   gap <- read.table(
#     gzfile(gaps_file),
#     header = FALSE,
#     stringsAsFactors = FALSE
#   )
#   
#   colnames(gap) <- c(
#     "bin",
#     "chr",
#     "start",
#     "end",
#     "ix",
#     "n",
#     "size",
#     "type",
#     "bridge"
#   )
#   
#   gaps <- GRanges(
#     seqnames = gap$chr,
#     ranges = IRanges(
#       start = gap$start + 1,  # UCSC BED-style 0-based -> GRanges 1-based
#       end = gap$end
#     ),
#     type = gap$type
#   )
#   assembly_gaps <- gaps[gaps$type %in% c(
#     "centromere",
#     "telomere",
#     "heterochromatin",
#     "short_arm",
#     "contig"
#   )]
#   
#   
#   gr_background <- setdiff(
#     gr_outside_peaks,
#     c(blacklist_gr, centromeres, assembly_gaps)
#   )
#   standard_chromosomes <- paste0('chr',seq_along(1:22))
#   bg <- as.data.frame(gr_background) %>% 
#     filter(seqnames%in%standard_chromosomes) %>% 
#     filter(width>=filter_small_than)
#   colnames(bg)<-c('bg_chr','bg_start','bg_end','bg_width','bg_strand')
#   return(bg)
# }



merge_intervals <- function(df){
  
  if(nrow(df) == 0)
    return(df)
  
  df <- df[order(df$start, df$end), ]
  
  starts <- c()
  ends <- c()
  
  cur_start <- df$start[1]
  cur_end <- df$end[1]
  
  if(nrow(df) > 1){
    for(i in 2:nrow(df)){
      
      if(df$start[i] <= cur_end + 1){
        
        cur_end <- max(cur_end, df$end[i])
        
      } else {
        
        starts <- c(starts, cur_start)
        ends <- c(ends, cur_end)
        
        cur_start <- df$start[i]
        cur_end <- df$end[i]
      }
    }
  }
  
  starts <- c(starts, cur_start)
  ends <- c(ends, cur_end)
  
  data.frame(
    start = starts,
    end = ends
  )
}


complement_intervals <- function(blocked, chr_length){
  
  bg_start <- c()
  bg_end <- c()
  
  prev <- 1
  
  if(nrow(blocked) > 0){
    
    for(i in seq_len(nrow(blocked))){
      
      if(prev < blocked$start[i]){
        
        bg_start <- c(bg_start, prev)
        bg_end <- c(bg_end, blocked$start[i]-1)
        
      }
      
      prev <- max(prev, blocked$end[i]+1)
    }
  }
  
  if(prev <= chr_length){
    
    bg_start <- c(bg_start, prev)
    bg_end <- c(bg_end, chr_length)
    
  }
  
  data.frame(
    start = bg_start,
    end = bg_end
  )
}

# simulate_background_fragments <- function(
#     background_regions,
#     lambda_per_kb = 0.01,
#     frag_len_out_peak_dens
# ) {
#   
#   # Convert to dataframe
#   standard_chromosomes <- paste0('chr',seq_along(1:22))
#   bg <- background_regions %>% 
#     filter(bg_chr%in%standard_chromosomes)
#   
#   
#   
#   # total background size per chromosome
#   chr_sizes <- bg %>%
#     group_by(bg_chr) %>%
#     summarise(total_bp = sum(bg_width))
#   
#   # 1. simulate fragment counts per chromosome
#   chr_sizes <- chr_sizes %>%
#     dplyr::mutate(
#       lambda = lambda_per_kb * total_bp / 10000,
#       n_frag = rpois(n(), lambda)
#     )
#   
#   fragments <- list()
#   
#   for (i in seq_len(nrow(chr_sizes))) {
#     
#     chr <- chr_sizes$bg_chr[i]
#     n <- chr_sizes$n_frag[i]
#     
#     if (n == 0)
#       next
#     
#     # background intervals for this chromosome
#     chr_bg <- bg %>%
#       dplyr::filter(bg_chr == chr) %>% 
#       dplyr::mutate(bg_start_100bb=bg_start-100,
#              bg_end_100bb=bg_end+100) %>% 
#       dplyr::mutate(bg_width_flank=bg_end_100bb-bg_start_100bb) %>% 
#       dplyr::filter(bg_width_flank>0)
#     
#     # probability proportional to interval width
#     idx <- sample(
#       seq_len(nrow(chr_bg)),
#       size = n,
#       replace = TRUE,
#       prob = chr_bg$bg_width_flank
#     )
#     
#     # sample a position within each selected interval
#     starts <- chr_bg$bg_start_100bb[idx] + sample.int(
#       max(chr_bg$bg_width_flank),
#       n,
#       replace = TRUE
#     ) %% chr_bg$bg_width_flank[idx]
#     
#     # 3. sample fragment 
#     frag_sampler <- function(my_frag_leng_dist,n) {
#       sample(
#         my_frag_leng_dist$x,
#         size = n,
#         replace = TRUE,
#         prob = my_frag_leng_dist$y
#       )
#     }
#     frag_size <- frag_sampler(my_frag_leng_dist=frag_len_out_peak_dens,n)
#     frag_size <- as.integer(round(frag_size))
#     ends <- starts + frag_size - 1
#     
#     
#     # remove fragments crossing background interval boundary
#     # valid <- ends <= chr_bg$bg_end
#     
#     fragments[[i]] <- data.frame(
#       fragment = paste0("bg_frag_", seq_len(n)),
#       fragment_chr = str_remove(chr,pattern = "chr"),
#       fragment_start = starts,
#       fragment_end = ends,
#       fragment_size = frag_size
#     )
#   }
#   
#   bind_rows(fragments)
# }


simulate_background_fragments <- function(
    background_regions,
    mean_pct_fragments_out,
    tot_fragments_in_peak,
    frag_len_out_peak_dens
) {
  
  
  pct_out_peak_fragments = rnorm(n = 1,mean = mean_pct_fragments_out,sd = 0.05)
  tot_frags_out <- round(tot_fragments_in_peak *  pct_out_peak_fragments / (1-pct_out_peak_fragments),0)
  
  
  
  # Convert to dataframe
  standard_chromosomes <- paste0('chr',seq_along(1:22))
  bg <- background_regions %>% 
    dplyr::filter(bg_chr%in%standard_chromosomes) %>% 
    dplyr::mutate(bg_start_offset=round(bg_start-100),0) %>% 
    dplyr::mutate(bg_end_offset=round(bg_end+100)) %>% 
    dplyr::mutate(bg_withd_offset=bg_end_offset-bg_start_offset)
  
  bg_bed <- data.table(bg)
  
  # Randomly sample bed file rows, proportional to the length of each range
  simulated.sites <- bg_bed[sample(.N, size=tot_frags_out, replace=TRUE, prob=bg_bed$bg_withd_offset)]
  
  # Randomly sample uniformly within each chosen range
  simulated.sites[, fragment_start := sample(bg_start_offset:bg_end_offset, size=1), by=1:dim(simulated.sites)[1]]
  bg=as.data.frame(simulated.sites)
  frag_sampler <- function(my_frag_leng_dist,n) {
    sample(
      my_frag_leng_dist$x,
      size = n,
      replace = TRUE,
      prob = my_frag_leng_dist$y
    ) %>% round(0)
  }
  
  
  bg <- bg %>% 
    dplyr::mutate(fragment_size=sample(frag_len_out_peak_dens,size = nrow(simulated.sites),replace = TRUE)) %>% 
    # dplyr::mutate(fragment_size=frag_sampler(my_frag_leng_dist=frag_len_out_peak_dens,n=nrow(simulated.sites))) %>% 
    dplyr::mutate(fragment_end=fragment_start+fragment_size-1) %>% 
    dplyr::mutate(fragment=paste0("bg_frag_", seq_len(nrow(simulated.sites)))) %>% 
    dplyr::rename(fragment_chr=bg_chr) %>% 
    dplyr::select(starts_with("fragment"))
  bg$fragment_chr <- gsub(pattern = 'chr',replacement = '',x=bg$fragment_chr )
  return(bg)
}



get_background_regions <- function(
    peak_fragments_df,
    chrom_sizes_file,
    blacklist_file,
    gaps_file,
    centromeres_file,
    filter_small_than = 150
){
  

  
  peak_region_df <- peak_fragments_df %>%
    dplyr::group_by(peak) %>%
    dplyr::summarise(
      start=min(fragment_start),
      end=max(fragment_end),
      chr=paste0("chr", unique(fragment_chr)),
      .groups="drop"
    )

  
  chrom_sizes <- chrom_sizes_file %>% 
    dplyr::rename(chr=V1,
                  length=V3) %>% 
    dplyr::mutate(chr=paste0('chr',chr))

  
  blacklist <-read.table(blacklist_file,col.names = c("chr","start","end"))
  
  
  ## BED -> 1-based
  blacklist$start <- blacklist$start + 1
  

  
  cytoband <- read.table(
    gzfile(centromeres_file),
    header=FALSE,
    sep="\t",
    stringsAsFactors=FALSE
  )
  
  colnames(cytoband) <- c(
    "chr",
    "start",
    "end",
    "name",
    "stain"
  )
  
  centromeres <- cytoband %>%
    filter(stain=="acen") %>%
    transmute(
      chr,
      start=start+1,
      end=end
    )

  
  gap <- read.table(
    gzfile(gaps_file),
    header=FALSE,
    stringsAsFactors=FALSE
  )
  
  colnames(gap) <- c(
    "bin",
    "chr",
    "start",
    "end",
    "ix",
    "n",
    "size",
    "type",
    "bridge"
  )
  
  assembly_gaps <- gap %>%
    filter(type %in% c(
      "centromere",
      "telomere",
      "heterochromatin",
      "short_arm",
      "contig"
    )) %>%
    transmute(
      chr,
      start=start+1,
      end=end
    )
  

  
  blocked <- bind_rows(
    peak_region_df,
    blacklist,
    centromeres,
    assembly_gaps
  )

  
  background_list <- list()
  
  for(chr in chrom_sizes$chr){
    
    chr_length <- chrom_sizes$length[chrom_sizes$chr==chr]
    
    chr_blocked <- blocked %>%
      dplyr::filter(chr==!!chr) %>%
      dplyr::select(start,end)
    
    if(nrow(chr_blocked)>0){
      
      chr_blocked <- merge_intervals(chr_blocked)
      
    }
    
    bg <- complement_intervals(
      chr_blocked,
      chr_length
    )
    
    if(nrow(bg)>0){
      
      bg$chr <- chr
      
      background_list[[chr]] <- bg
      
    }
    
  }
  
  bg <- bind_rows(background_list) %>%
    dplyr::mutate(
      width=end-start+1
    ) %>%
    dplyr::filter(
      chr %in% paste0("chr",1:22),
      width>=filter_small_than
    ) %>%
    dplyr::select(
      bg_chr=chr,
      bg_start=start,
      bg_end=end,
      bg_width=width
    )
  
  bg$bg_strand <- "*"
  bg
}

get_epigenetic_activity<- function(activity,epistate){
  programs <- activity[[epistate]]
  return(programs)
}



add_sparsity_all<- function(real_df,dropout_rate=0.7){
  N_obs <- nrow(real_df)
  N_real_zeros <- real_df %>%
    filter(status == 0) %>%
    nrow()
  N_real_ones <- real_df %>%
    filter(status == 1) %>%
    nrow()
  # current real zero proportion
  real_zero_pct <- N_real_zeros / N_obs
  # additional dropout needed among the existing ones
  remaining_dropout <- dropout_rate - real_zero_pct
  
  # number of accessible peaks to remove
  n_drop <- round(remaining_dropout * N_real_ones)
  
  # avoid negative values if already above target sparsity
  if (n_drop < 0) {
    n_drop <- 0
  }
  
  # sample ONLY currently accessible peaks
  set.seed(123)
  drop_idx <- sample(
    which(real_df$status == 1),
    size = n_drop,
    replace = FALSE
  )
  
  # copy dataframe
  df_dropout <- real_df
  
  # convert dropout peaks to zeros
  df_dropout$status[drop_idx] <- 0
  return(df_dropout)
}

rtrunc_weibull <- function(n, shape=84, scale=0.85, lower = 0, upper = 0.97) {
  x <- numeric(n)
  i <- 1
  
  while (i <= n) {
    y <- rweibull(n - i + 1, shape = shape, scale = scale)
    y <- y[y >= lower & y <= upper]
    if (length(y) > 0) {
      k <- min(length(y), n - i + 1)
      x[i:(i + k - 1)] <- y[1:k]
      i <- i + k
    }
  }
  
  x
}

add_sparsity_per_cell <- function(real_df,weibull_scale=0.9){
  
  cells = unique(real_df$cell_id)
  cell_sparsity <- tibble(
    cell_id = cells,
    sparsity = rtrunc_weibull(length(cells), scale = weibull_scale)
  )
  
  real_df_after <- real_df %>%
    dplyr::left_join(cell_sparsity, by = "cell_id") %>%
    dplyr::mutate(
      status = if_else(
        status == 1,
        rbinom(n(), 1, 1 - sparsity),
        status
      )
    ) %>%
    dplyr::select(-sparsity)
  return(real_df_after)
}

library(dplyr)
library(purrr)
library(tibble)

convert_activity_list<-function(activity_list){
  activity_df <- imap_dfr(activity_list, function(cell, cell_name) {
    
    imap_dfr(cell, function(sign_values, sign_name) {
      
      tibble(
        mutant = cell_name,
        pathway = sign_name,
        # pathway = names(sign_values),
        activity = as.numeric(sign_values)
      )
      
    })
  })
  
  return(activity_df)
}

sample_fragments_for_peak_vec_allele <- function(
    peak_id,
    peak_chr,
    peak_from,
    peak_to,
    flank = 100,
    lambda = 0.2,
    fragment_len_dist = NULL,
    available_alleles,
    max_fragments = 50,
    cell_id
) {

  n <- length(peak_id)

  stopifnot(
    length(peak_from) == n,
    length(peak_to) == n
  )

  peak_size <- peak_to - peak_from
  max_total <- peak_size + 2 * flank

  # vectorized fragment count sampling
  n_frag_vec <- rztpois(n = n,lambda = lambda)

  n_frag_vec <- pmin(n_frag_vec, max_fragments)

  results <- vector("list", n)

  for (k in seq_len(n)) {
    n_frag <- n_frag_vec[k]

    if (n_frag == 0 || length(available_alleles) == 0) {
      next
    }

    allele_results <- list()

    for (allele in as.character(available_alleles)) {

      # sample all fragment sizes simultaneously
      if (!is.null(fragment_len_dist)) {

        sizes <- sample(
          fragment_len_dist,
          size = n_frag,
          replace = TRUE
        )

      } else {

        sizes <- sample_fragment_size(n_frag)

      }

      # sample all gaps simultaneously
      gaps <- numeric(n_frag)

      if (n_frag > 1) {
        gaps[2:n_frag] <- sample_gap_size(n_frag - 1)
      }

      # cumulative occupied space
      occupied <- cumsum(sizes + gaps)

      keep <- occupied <= max_total[k]

      if (!any(keep)) {
        next
      }

      sizes <- sizes[keep]
      gaps  <- gaps[keep]
      print(peak_id[k])
      print(gaps)
      print(sizes)
      # vectorized fragment placement
      peaks_frags_df <- place_fragments_in_peak(fragments_sizes = sizes,
                                                gaps_sizes = gaps,
                                                peak_id = peak_id[k],peak_from = peak_from[k],
                                                peak_to = peak_to[k],sd_peak_center=100)

      allele_results[[allele]] <- peaks_frags_df %>% mutate(fragment_allele=allele) %>%
        mutate(peak_id = peak_id[k],
               peak_chr=peak_chr,
               peak_from = peak_from[k],
               peak_to = peak_to[k],
               fragment_size = sizes,
               cell_id=cell_id,
               fragment_chr=peak_chr
        ) %>%
        dplyr::rename(fragment_start=from,
                      fragment_end=to)

    }

    results[[k]] <- do.call("rbind",allele_results)
    # results_pre_sparsity <- do.call("rbind",allele_results)
    
  }

  do.call("rbind",results)
}

# library(dplyr)
# library(GenomeInfoDb)
# library(BSgenome.Hsapiens.UCSC.hg38)
# 
# # chromosome lengths
# realtive_to_absolute <- function(df,genome){
#   if (genome=='hg38'){
#     seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
#     
#     # cumulative offsets
#     offsets <- cumsum(c(0, seqlens[names(seqlens) != "chrY"][-length(seqlens[names(seqlens) != "chrY"])]))
#     
#     chr_offsets <- data.frame(
#       chrom = names(offsets),
#       offset = offsets
#     )
#     
#     # your data
#     seg <- df
#     
#     # convert
#     seg_abs <- seg %>%
#       left_join(chr_offsets, by = "chrom") %>%
#       mutate(
#         abs_start = start + offset,
#         abs_end   = end + offset
#       )
#   }
#   return(seg_abs)
# }


evaluate_peaks <- function(truth, called, min_overlap = 1) {
  truth_gr <- GRanges(
    seqnames = truth$chrom,
    ranges = IRanges(start = truth$start, end = truth$end))
  
  called_gr <- GRanges(
    seqnames = called$chrom,
    ranges = IRanges(start = called$start, end = called$end)
  )
  
  hits <- findOverlaps(
    called_gr,
    truth_gr,
    minoverlap = min_overlap
  )
  colnames(truth)<- paste0('simulated_',colnames(truth))
  colnames(called)<- paste0('inferred_',colnames(called))
  
  matched <- data.frame(
    truth_idx  = subjectHits(hits),
    called_idx = queryHits(hits)
  )
  
  matched <- cbind(
    truth[matched$truth_idx, ],
    called[matched$called_idx, ]
  )
  
  TP_called <- n_distinct(matched$inferred_name)
  TP_truth  <- n_distinct(matched$simulated_peak_id)
  
  # Peak-level TP count
  TP <- TP_called
  
  # Total peaks
  N_called <- nrow(called)
  N_truth  <- nrow(truth)
  
  
  # False positives and false negatives
  FP <- N_called - TP_called
  FN <- N_truth - TP_truth
  
  # Metrics
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  
  f1 <- ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
  
  peak_metrics <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
  
  return(list('peak_matching_df'=matched,'metrics'=peak_metrics))
 
}


evaluate_peaks_new <- function(truth, called, min_overlap = 10) {
  
  truth_gr <- GRanges(
    seqnames = truth$peak_chr,
    ranges = IRanges(start = truth$peak_from, end = truth$peak_to),
    pct_cells = truth$pct_cells,
  )
  
  called_gr <- GRanges(
    seqnames = called$chrom,
    ranges = IRanges(start = called$start, end = called$end)
  )
  
  hits <- findOverlaps(
    called_gr,
    truth_gr,
    minoverlap = min_overlap
  )
  
  colnames(truth)  <- paste0("simulated_", colnames(truth))
  colnames(called) <- paste0("inferred_", colnames(called))

  matched_idx <- data.frame(
    truth_idx  = subjectHits(hits),
    called_idx = queryHits(hits)
  )
  
  tp_df <- cbind(
    truth[matched_idx$truth_idx, , drop = FALSE],
    called[matched_idx$called_idx, , drop = FALSE]
  )
  tp_df$status <- "TP"

  fn_idx <- setdiff(seq_len(nrow(truth)), unique(subjectHits(hits)))
  
  fn_df <- truth[fn_idx, , drop = FALSE]
  fn_df$status <- "FN"

  fp_idx <- setdiff(seq_len(nrow(called)), unique(queryHits(hits)))
  
  fp_df <- called[fp_idx, , drop = FALSE]
  fp_df$status <- "FP"
  
  ## Make dataframes compatible for rbind()
  all_cols <- union(colnames(tp_df), union(colnames(fn_df), colnames(fp_df)))
  
  add_missing <- function(df, cols) {
    miss <- setdiff(cols, colnames(df))
    df[miss] <- NA
    df[, cols, drop = FALSE]
  }
  
  tp_df <- add_missing(tp_df, all_cols)
  fn_df <- add_missing(fn_df, all_cols)
  fp_df <- add_missing(fp_df, all_cols)
  
  final_df <- rbind(tp_df, fn_df, fp_df)
  
  ## Metrics
  TP_called <- dplyr::n_distinct(tp_df$inferred_name)
  TP_truth  <- dplyr::n_distinct(tp_df$simulated_peak)
  
  TP <- TP_called
  
  N_called <- nrow(called)
  N_truth  <- nrow(truth)
  
  FP <- N_called - TP_called
  FN <- N_truth - TP_truth
  
  precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
  recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
  
  f1 <- ifelse(
    precision + recall == 0,
    0,
    2 * precision * recall / (precision + recall)
  )
  
  peak_metrics <- data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
  
  return(list(
    peak_matching_df = tp_df,
    peak_classification_df = final_df,
    metrics = peak_metrics,
    min_overlap = min_overlap
  ))
}


evaluate_peaks_extended <- function(truth, called, min_overlap = 10) {
  
  # Harmonise chromosome names
  truth$peak_chr <- ifelse(
    grepl("^chr", truth$peak_chr),
    truth$peak_chr,
    paste0("chr", truth$peak_chr)
  )
  
  called$chrom <- ifelse(
    grepl("^chr", called$chrom),
    called$chrom,
    paste0("chr", called$chrom)
  )
  
  # Remove duplicate ground-truth regions
  truth <- truth %>%
    distinct(peak_chr, peak_from, peak_to, .keep_all = TRUE)
  
  # BED/narrowPeak: 0-based, half-open coordinates.
  # Convert to GRanges: 1-based, closed coordinates.
  # Adjust the truth conversion if its coordinates are
  # already 1-based and inclusive.
  truth_gr <- GRanges(
    seqnames = truth$peak_chr,
    ranges = IRanges(
      start = truth$peak_from + 1L,
      end = truth$peak_to
    )
  )
  
  called_gr <- GRanges(
    seqnames = called$chrom,
    ranges = IRanges(
      start = called$start + 1L,
      end = called$end
    )
  )
  
  hits <- findOverlaps(
    called_gr,
    truth_gr,
    minoverlap = min_overlap
  )
  
  called_matched <- unique(queryHits(hits))
  truth_matched <- unique(subjectHits(hits))
  # Map simulated regulatory regions to MACS3-called peaks
  matched_idx <- data.frame(
    truth_idx = subjectHits(hits),
    called_idx = queryHits(hits)
  )
  
  # Calculate the number of overlapping base pairs
  matched_idx$overlap_bp <- width(
    pintersect(
      truth_gr[matched_idx$truth_idx],
      called_gr[matched_idx$called_idx]
    )
  )
  
  # Create the mapping dataframe
  peak_mapping_df <- cbind(
    truth[matched_idx$truth_idx, , drop = FALSE] %>%
      ungroup() %>%
      rename_with(~ paste0("simulated_", .x)),
    
    called[matched_idx$called_idx, , drop = FALSE] %>%
      ungroup() %>%
      rename_with(~ paste0("inferred_", .x))
  ) %>%
    mutate(
      overlap_bp = matched_idx$overlap_bp
    ) %>%
    as.data.frame()
  
  peak_mapping_df <- peak_mapping_df %>%
    mutate(
      simulated_width = simulated_peak_to - simulated_peak_from,
      inferred_width = inferred_end - inferred_start,
      
      fraction_simulated_covered =
        overlap_bp / simulated_width,
      
      fraction_inferred_covered =
        overlap_bp / inferred_width
    ) %>% 
    group_by(inferred_name) %>%
    mutate(
      n_simulated_regions = n_distinct(simulated_peak)
    ) %>%
    ungroup()
  
  
  n_called <- length(called_gr)
  n_truth <- length(truth_gr)
  
  n_supported <- length(called_matched)
  n_recovered <- length(truth_matched)
  
  precision <- if (n_called > 0) {
    n_supported / n_called
  } else NA_real_
  
  recall <- if (n_truth > 0) {
    n_recovered / n_truth
  } else NA_real_
  
  f1 <- if (is.na(precision) || is.na(recall)) {
    NA_real_
  } else if (precision + recall > 0) {
    2 * precision * recall / (precision + recall)
  } else 0
  
  # Base-pair-level overlap
  truth_union <- reduce(truth_gr)
  called_union <- reduce(called_gr)
  
  intersection_bp <- sum(
    width(intersect(truth_union, called_union))
  )
  
  union_bp <- sum(
    width(union(truth_union, called_union))
  )
  
  jaccard <- if (union_bp > 0) {
    intersection_bp / union_bp
  } else NA_real_
  
  # Number of overlapping regions per peak
  called_overlap_counts <- tabulate(
    queryHits(hits), nbins = n_called
  )
  
  # Number of overlapping peaks per truth region
  truth_overlap_counts <- tabulate(
    subjectHits(hits), nbins = n_truth
  )
  
  called_results <- called %>%
    mutate(
      recovered_truth_regions = called_overlap_counts,
      supported = recovered_truth_regions > 0,
      merges_regions = recovered_truth_regions > 1
    )
  
  truth_results <- truth %>%
    mutate(
      overlapping_called_peaks = truth_overlap_counts,
      recovered = overlapping_called_peaks > 0,
      split_across_peaks = overlapping_called_peaks > 1
    )
  
  metrics <- data.frame(
    n_called = n_called,
    n_truth = n_truth,
    n_supported_called = n_supported,
    n_recovered_truth = n_recovered,
    precision = precision,
    recall = recall,
    f1 = f1,
    jaccard = jaccard,
    n_peaks_merging_regions =
      sum(called_overlap_counts > 1),
    n_regions_overlapping_multiple_peaks =
      sum(truth_overlap_counts > 1)
  )
  
  list(
    metrics = metrics,
    peak_mapping_df = peak_mapping_df,
    truth_results = truth_results,
    called_results = called_results,
    hits = hits
  )
}

library(ggplot2)
library(dplyr)
library(patchwork)

plot_peak_example <- function(
    peak_name,
    results,
    fragments,
    flank = 500
) {
  
  called <- results$called_results %>% ungroup()
  truth <- results$truth_results %>% ungroup()
  
  # Select MACS3 peak
  selected <- called %>%
    filter(name == peak_name)
  
  if (nrow(selected) != 1) {
    stop("Peak name not found or not unique.")
  }
  
  chr <- selected$chrom[1]
  
  region_start <- max(0L, selected$start[1] - flank)
  region_end <- selected$end[1] + flank
  
  # Simulated regulatory regions
  truth_local <- truth %>%
    filter(
      peak_chr == chr,
      peak_from < region_end,
      peak_to > region_start
    )
  
  # All MACS3 peaks within the selected interval
  called_local <- called %>%
    filter(
      chrom == chr,
      start < region_end,
      end > region_start
    )
  
  # Sequenced fragments overlapping the region
  frags_local <- fragments %>%
    filter(
      fragment_chr == chr,
      fragment_start < region_end,
      fragment_end > region_start
    )
  
  # Calculate coverage only for the selected interval
  window_width <- as.integer(region_end - region_start)
  
  coverage_values <- numeric(window_width)
  
  if (nrow(frags_local) > 0) {
    
    # Clip fragments to the plotting window
    frag_start <- pmax(
      as.integer(frags_local$fragment_start),
      region_start
    )
    
    frag_end <- pmin(
      as.integer(frags_local$fragment_end),
      region_end
    )
    
    valid <- frag_end > frag_start
    
    if (any(valid)) {
      ir <- IRanges(
        start = frag_start[valid] - region_start + 1L,
        end = frag_end[valid] - region_start
      )
      
      coverage_values <- as.numeric(
        IRanges::coverage(ir, width = window_width)
      )
    }
  }
  
  coverage_df <- data.frame(
    position = seq.int(
      region_start,
      region_end - 1L
    ),
    coverage = coverage_values
  )
  
  # Track 1: simulated regulatory regions
  p_truth <- ggplot(truth_local) +
    geom_segment(
      aes(
        x = peak_from,
        xend = peak_to,
        y = 1,
        yend = 1
      ),
      linewidth = 5,
      colour = "#3A9D7A"
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end)
    ) +
    labs(y = "Simulated") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  # Track 2: MACS3 peaks
  p_called <- ggplot(called_local) +
    geom_segment(
      aes(
        x = start,
        xend = end,
        y = 1,
        yend = 1
      ),
      linewidth = 6,
      colour = "#4472C4"
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end)
    ) +
    labs(y = "MACS3") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank()
    )
  
  # Track 3: aggregated fragment coverage
  p_coverage <- ggplot(
    coverage_df,
    aes(x = position, y = coverage)
  ) +
    geom_area(
      fill = "#8B5CF6",
      alpha = 0.65
    ) +
    geom_line(
      colour = "#6D28D9",
      linewidth = 0.3
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end)
    ) +
    labs(
      x = paste0(chr, " genomic position (bp)"),
      y = "Fragment coverage"
    ) +
    theme_minimal()
  
  # Combine tracks
  p_truth / p_called / p_coverage +
    plot_layout(heights = c(1, 1, 3)) +
    plot_annotation(
      title = paste("Genomic example:", peak_name)
    )
}

library(dplyr)
library(ggplot2)
library(patchwork)
library(IRanges)




my_ggplot_theme <- function () 
{
  theme_light(base_size=10) +
    theme(legend.key.size=unit(0.3, "cm"),
          panel.background=element_rect(fill="white"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.title=element_text(size=10),
          text=element_text(size=10),
          plot.title=element_text(size=12))
}

plot_genomic_region <- function(
    chr,
    region_start,
    region_end,
    results,
    fragments
) {
  
  stopifnot(
    region_end > region_start,
    region_start >= 0
  )
  
  # BED-style coordinates: 0-based, half-open
  chr <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  
  called <- results$called_results %>% ungroup()
  truth <- results$truth_results %>% ungroup()
  
  # 1. Select simulated regulatory regions
  truth_local <- truth %>%
    filter(
      peak_chr == chr,
      peak_from < region_end,
      peak_to > region_start
    )
  
  # 2. Select all MACS3 peaks in the region
  called_local <- called %>%
    filter(
      chrom == chr,
      start < region_end,
      end > region_start
    ) %>%
    mutate(
      midpoint = (
        pmax(start, region_start) +
          pmin(end, region_end)
      ) / 2
    )
  
  # 3. Select sequenced fragments in the region
  frags_local <- fragments %>%
    filter(
      fragment_chr == chr,
      fragment_start < region_end,
      fragment_end > region_start,
      !is.na(fragment_type)
    )
  
  # 4. Calculate coverage by fragment type
  window_width <- as.integer(region_end - region_start)
  
  fragment_types <- unique(
    as.character(frags_local$fragment_type)
  )
  
  coverage_list <- lapply(fragment_types, function(type) {
    
    frags_type <- frags_local %>%
      filter(fragment_type == type)
    
    frag_start <- pmax(
      as.integer(frags_type$fragment_start),
      region_start
    )
    
    frag_end <- pmin(
      as.integer(frags_type$fragment_end),
      region_end
    )
    
    valid <- frag_end > frag_start
    coverage_values <- numeric(window_width)
    
    if (any(valid)) {
      
      ir <- IRanges::IRanges(
        start = frag_start[valid] - region_start + 1L,
        end = frag_end[valid] - region_start
      )
      
      coverage_values <- as.numeric(
        IRanges::coverage(
          ir,
          width = window_width
        )
      )
    }
    
    data.frame(
      position = seq.int(region_start, region_end - 1L),
      coverage = coverage_values,
      fragment_type = type
    )
  })
  
  coverage_df <- if (length(coverage_list) > 0) {
    bind_rows(coverage_list)
  } else {
    data.frame(
      position = seq.int(region_start, region_end - 1L),
      coverage = 0,
      fragment_type = "No fragments"
    )
  }
  
  # Consistent fragment type colours
  fragment_colours <- c(
    peak = "#8B5CF6",
    background = "#BDBDBD"
  )
  
  # Assign colours to any additional fragment types
  extra_types <- setdiff(
    unique(coverage_df$fragment_type),
    names(fragment_colours)
  )
  
  if (length(extra_types) > 0) {
    extra_colours <- setNames(
      grDevices::hcl.colors(
        length(extra_types),
        palette = "Dark 3"
      ),
      extra_types
    )
    
    fragment_colours <- c(
      fragment_colours,
      extra_colours
    )
  }
  
  
  # 6. Simulated regulatory regions
  p_truth <- ggplot(truth_local) +
    geom_segment(
      aes(
        x = peak_from,
        xend = peak_to,
        y = 1,
        yend = 1
      ),
      linewidth = 5,
      colour = "#3A9D7A"
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      ylim = c(0.5, 1.5),
      expand = FALSE
    ) +
    labs(y = "Simulated regions",x="") +
    my_ggplot_theme
  
  # 7. MACS3 peaks with labels
  p_called <- ggplot(called_local) +
    geom_segment(
      aes(
        x = start,
        xend = end,
        y = 1,
        yend = 1
      ),
      linewidth = 6,
      colour = "#4472C4"
    ) +
    geom_text(
      aes(
        x = midpoint,
        y = 1.35,
        label = name
      ),
      size = 3,
      colour = "#1E3A8A",
      check_overlap = TRUE
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      ylim = c(0.5, 1.7),
      expand = FALSE
    ) +
    labs(y = "MACS3 peaks",x="") +
    my_ggplot_theme
  
  # 8. Stacked fragment coverage
  p_coverage <- ggplot(
    coverage_df,
    aes(
      x = position,
      y = coverage,
      fill = fragment_type
    )
  ) +
    geom_area(
      position = "stack",
      alpha = 0.85,
      linewidth = 0
    ) +
    scale_fill_manual(
      values = fragment_colours
    ) +
    coord_cartesian(
      xlim = c(region_start, region_end),
      expand = FALSE
    ) +
    facet_wrap(~fragment_type,ncol=1)+
    scale_x_continuous(
      labels = scales::label_comma()
    ) +
    labs(
      x = paste0(chr, " genomic position (bp)"),
      y = "Fragment coverage",
      fill = "Fragment type"
    ) +
    my_ggplot_theme +
    theme(
      legend.position = "bottom"
    )
  
  # 9. Combine tracks
  p_truth / p_called / p_coverage +
    plot_layout(
      heights = c(1, 1.5, 3)
    ) +
    plot_annotation(
      title = paste0(
        chr, ":",
        format(region_start, big.mark = ","),
        "-",
        format(region_end, big.mark = ",")
      )
    )
}
