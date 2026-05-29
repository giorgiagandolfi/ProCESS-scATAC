library(extraDistr)
library(data.table)
library(Matrix)
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
sample_fragment_size <- function(n) {
  
  comp <- sample(
    x = 1:4,
    size = n, # change to n fragments
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

place_fragments <- function(start, fragments_sizes,gaps_sizes) {
  n <- length(fragments_sizes)
  from <- numeric(n)
  to <- numeric(n)
  
  current <- start
  
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
      
      # vectorized fragment placement
      peaks_frags_df <- place_fragments_in_peak(fragments_sizes = sizes,
                                                gaps_sizes = gaps,
                                                peak_id = peak_id[k],peak_from = peak_from[k],
                                                peak_to = peak_to[k])
      
      allele_results[[allele]] <- peaks_frags_df %>% mutate(allele=paste0('allele_',allele)) %>% 
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
    
    results[[k]] <- data.table::rbindlist(
      allele_results,
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  data.table::rbindlist(
    results,
    use.names = TRUE,
    fill = TRUE
  )
}


