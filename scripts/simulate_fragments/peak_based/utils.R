# library(extraDistr)
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
    
    results[[k]] <- do.call("rbind",allele_results)
  }
  
  do.call("rbind",results)
}


get_background_regions <- function(peak_df,genome,gaps_file,centromeres_file){
  #### peak_df is a dataframe with all the peaks for that cell
  #### like this one
  # peak_id peak_chr peak_start  peak_end
  # 1        chr1-1006219-1006719     chr1    1006219   1006719
  # 2        chr1-1068981-1069481     chr1    1068981   1069481
  # 3      chr1-10694533-10695033     chr1   10694533  10695033
  # 4    chr1-108559703-108560203     chr1  108559703 108560203
  # 5    chr1-109763665-109764165     chr1  109763665 109764165
  # 6    chr1-111204245-111204745     chr1  111204245 111204745
  gr_peak_union = GRanges(
    seqnames = peak_df$chr,
    ranges = IRanges(start = peak_df$from, end = peak_df$to),
  )
  if (genome=='hg38'){
    gr_genome <- GRanges(
      seqnames = names(seqlengths(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)),
      ranges = IRanges(
        start = 1,
        end = seqlengths(BSgenome.Hsapiens.UCSC.hg38)
      ))
  } else if (genome=='hg19'){
    gr_genome <- GRanges(
      seqnames = names(seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)),
      ranges = IRanges(
        start = 1,
        end = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
      ))
  }
  gr_outside_peaks <- setdiff(gr_genome, gr_peak_union)
  
  blacklist_gr <- BiocIO::import(
    "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
  )
  
  blacklist_gr <- keepStandardChromosomes(
    blacklist_gr,
    pruning.mode = "coarse"
  )
  
  cytoband <- read.table(
    gzfile(centromeres_file), ## must be a gz
    header = FALSE,
    sep = '\t',
    col.names = c("chr", "start", "end", "name", "stain")
  )
  centromeres <- GRanges(
    seqnames = cytoband$chr[cytoband$stain == "acen"],
    ranges = IRanges(
      start = cytoband$start[cytoband$stain == "acen"] + 1, # BED is 0-based
      end = cytoband$end[cytoband$stain == "acen"]
    )
  )
  
  gap <- read.table(
    gzfile(gaps_file),
    header = FALSE,
    stringsAsFactors = FALSE
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
  
  gaps <- GRanges(
    seqnames = gap$chr,
    ranges = IRanges(
      start = gap$start + 1,  # UCSC BED-style 0-based -> GRanges 1-based
      end = gap$end
    ),
    type = gap$type
  )
  assembly_gaps <- gaps[gaps$type %in% c(
    "centromere",
    "telomere",
    "heterochromatin",
    "short_arm",
    "contig"
  )]
  
  
  gr_background <- setdiff(
    gr_outside_peaks,
    c(blacklist_gr, centromeres, assembly_gaps)
  )
  standard_chromosomes <- paste0('chr',seq_along(1:22))
  bg <- as.data.frame(gr_background) %>% 
    filter(seqnames%in%standard_chromosomes)
  colnames(bg)<-c('bg_chr','bg_start','bg_end','bg_width','bg_strand')
  return(bg)
}

simulate_background_fragments <- function(
    background_regions,
    lambda_per_kb = 0.01,
    frag_len_out_peak_dens
) {
  
  # Convert to dataframe
  standard_chromosomes <- paste0('chr',seq_along(1:22))
  bg <- background_regions %>% 
    filter(bg_chr%in%standard_chromosomes) %>% 
    filter(bg_width>=150)
  
  
  
  # total background size per chromosome
  chr_sizes <- bg %>%
    group_by(bg_chr) %>%
    summarise(total_bp = sum(bg_width))
  
  # 1. simulate fragment counts per chromosome
  chr_sizes <- chr_sizes %>%
    mutate(
      lambda = lambda_per_kb * total_bp / 10000,
      n_frag = rpois(n(), lambda)
    )
  
  fragments <- list()
  
  for (i in seq_len(nrow(chr_sizes))) {
    
    chr <- chr_sizes$bg_chr[i]
    n <- chr_sizes$n_frag[i]
    
    if (n == 0)
      next
    
    # background intervals for this chromosome
    chr_bg <- bg %>%
      filter(bg_chr == chr)
    
    # sample interval
    interval_idx <- sample(
      seq_len(nrow(chr_bg)),
      size = n,
      replace = TRUE
    )
    
    sampled_intervals <- chr_bg[interval_idx, ]
    # sampled_intervals <- distinct(sampled_intervals)
    # 2. sample start uniformly in background intervals
    starts <- mapply(
      function(s, e) {
        sample(seq(s, e), 1)
      },
      sampled_intervals$bg_start,
      sampled_intervals$bg_end
    )
    
    # 3. sample fragment 
    frag_sampler <- function(my_frag_leng_dist,n) {
      sample(
        my_frag_leng_dist$x,
        size = n,
        replace = TRUE,
        prob = my_frag_leng_dist$y
      )
    }
    frag_size <- frag_sampler(my_frag_leng_dist=frag_len_out_peak_dens,n)
    frag_size <- as.integer(round(frag_size))
    ends <- starts + frag_size - 1
    
    
    # remove fragments crossing background interval boundary
    valid <- ends <= sampled_intervals$bg_end
    
    fragments[[i]] <- data.frame(
      frag_id = paste0("bg_frag_", seq_len(n))[valid],
      frag_chr = chr,
      frag_start = starts[valid],
      frag_end = ends[valid],
      frag_size = frag_size[valid]
    )
  }
  
  bind_rows(fragments)
}
get_epigenetic_activity<- function(activity,mutant,clone){
  programs <- activity[[mutant]][[clone]]
  return(programs)
}
