# library(extraDistr)
library(data.table)
library(Matrix)
library(markovchain)
library(ComplexHeatmap)
library(grid)
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

simulate_background_fragments <- function(
    background_regions,
    lambda_per_kb = 0.01,
    frag_len_out_peak_dens
) {
  
  # Convert to dataframe
  standard_chromosomes <- paste0('chr',seq_along(1:22))
  bg <- background_regions %>% 
    filter(bg_chr%in%standard_chromosomes)
  
  
  
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
      filter(bg_chr == chr) %>% 
      mutate(bg_start_100bb=bg_start-100,
             bg_end_100bb=bg_end+100) %>% 
      mutate(bg_width_flank=bg_end_100bb-bg_start_100bb) %>% 
      filter(bg_width_flank>0)
    
    # probability proportional to interval width
    idx <- sample(
      seq_len(nrow(chr_bg)),
      size = n,
      replace = TRUE,
      prob = chr_bg$bg_width_flank
    )
    
    # sample a position within each selected interval
    starts <- chr_bg$bg_start_100bb[idx] + sample.int(
      max(chr_bg$bg_width_flank),
      n,
      replace = TRUE
    ) %% chr_bg$bg_width_flank[idx]
    
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
    # valid <- ends <= chr_bg$bg_end
    
    fragments[[i]] <- data.frame(
      fragment = paste0("bg_frag_", seq_len(n)),
      fragment_chr = str_remove(chr,pattern = "chr"),
      fragment_start = starts,
      fragment_end = ends,
      fragment_size = frag_size
    )
  }
  
  bind_rows(fragments)
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
    group_by(peak) %>%
    summarise(
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
      filter(chr==!!chr) %>%
      select(start,end)
    
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
    mutate(
      width=end-start+1
    ) %>%
    filter(
      chr %in% paste0("chr",1:22),
      width>=filter_small_than
    ) %>%
    select(
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

get_epigenetic_Q <- function(epigenetic_rates,epistate_name){
  states <- names(epigenetic_rates)
  
  # Initialize matrix
  mat <- matrix(0,
                nrow = length(states),
                ncol = length(states),
                dimnames = list(states, states))
  
  # Fill in switching rates
  for (from in states) {
    for (to in names(epigenetic_rates[[from]])) {
      if (!(to %in% c("duplication", "death"))) {
        mat[from, to] <- epigenetic_rates[[from]][[to]]
      }
    }
  }
  diag(mat) <- -rowSums(mat - diag(diag(mat)))
  plot_mat=Heatmap(mat,name = paste0("Transition matrix ",epistate_name),
                   show_heatmap_legend = FALSE,
                   col = colorRampPalette(c("white", "deeppink4"))(10),cluster_rows = F,cluster_columns = F,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
                   })
  plot_mat <- grid.grabExpr(
    draw(plot_mat)
  )
  ctmc <- new("ctmc", states = states, generator = mat, name = "epistates")
  ss<-steadyStates(ctmc)
  return(list("q_mat"=mat,"plot_q"=plot_mat,"steady_states"=ss))
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
        epistate = cell_name,
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

sample_shape <- function(nodes) {
  nodes %>%
    dplyr::mutate(shape_label = ifelse(is.na(.data$sample), "N/A",
                                       .data$sample)) %>%
    dplyr::pull(.data$shape_label)
}
plot_forest_with_bar <- function(
    forest,
    highlight_sample = NULL,
    color_map = NULL,
    alpha_function = NULL,
    shape_label_function = sample_shape,
    color_label_function = NULL,
    
    epistate_colors = c(
      "E1" = "forestgreen",
      "E2" = "goldenrod",
      "E3" = "orchid2"
    ),
    
    mutant_colors = c(
      "G1" = "coral2",
      "G2" = "turquoise4",
      "G3" = "darkorange"
    ),
    
    annotation_height = 0.35
) {
  
  
  if (!inherits(forest, "Rcpp_SampleForest") &&
      !inherits(forest, "Rcpp_PhylogeneticForest")) {
    
    stop('The parameter "forest" is not a ProCESS forest.')
  }
  
  
  
  forest_data <- forest$get_nodes()
  
  
  
  if (is.null(color_label_function)) {
    
    forest_data <- forest_data %>%
      ProCESS:::add_species_col("color_label")
    
    if (is.null(color_map)) {
      color_map <- get_species_colors(forest)
    }
    
  } else {
    
    if (is.null(color_map)) {
      stop(
        '"color_map" is mandatory when ',
        '"color_label_function" is specified.'
      )
    }
    
    forest_data[["color_label"]] <-
      color_label_function(forest_data)
  }
  
  
  forest_data$color_label <- factor(
    forest_data$color_label,
    levels = names(color_map)
  )
  
  
  
  if (!is.null(shape_label_function)) {
    
    forest_data[["shape_label"]] <-
      shape_label_function(forest_data)
    
  } else {
    
    forest_data[["shape_label"]] <- NA
  }
  
  
  
  if (is.null(alpha_function)) {
    
    forest_data <- forest_data %>%
      dplyr::mutate(
        alpha_level = 1
      )
    
  } else {
    
    forest_data <- forest_data %>%
      dplyr::mutate(
        alpha_level = alpha_function(.)
      )
  }
  
  
  
  if (nrow(forest_data) == 0) {
    
    warning("The forest does not contain any node")
    
    return(
      ggplot2::ggplot()
    )
  }
  
  
  
  forest_data <- forest_data %>%
    
    dplyr::as_tibble() %>%
    
    dplyr::rename(
      from = .data$ancestor,
      to   = .data$cell_id
    ) %>%
    
    dplyr::select(
      .data$from,
      .data$to,
      .data$sample,
      
      # Keep both annotations
      .data$epistate,
      .data$mutant,
      
      .data$shape_label,
      .data$color_label,
      .data$birth_time,
      .data$alpha_level
    )
  
  
  
  first_cell <- forest_data %>%
    dplyr::filter(
      .data$birth_time == 0
    )
  
  
  
  forest_data <- forest_data %>%
    
    dplyr::add_row(
      
      from = NA,
      to   = NA,
      
      color_label =
        first_cell[1, ]$color_label,
      
      epistate =
        first_cell[1, ]$epistate,
      
      mutant =
        first_cell[1, ]$mutant,
      
      shape_label =
        first_cell[1, ]$shape_label,
      
      birth_time = 0
    ) %>%
    
    dplyr::mutate(
      
      from = ifelse(
        is.na(.data$from),
        "WT",
        .data$from
      ),
      
      to = ifelse(
        is.na(.data$to),
        "WT",
        .data$to
      ),
      
      sample = ifelse(
        is.na(.data$sample),
        "N/A",
        as.character(.data$sample)
      ),
      
      highlight = FALSE
    )
  
  
  
  if (!is.null(highlight_sample)) {
    
    highlight <- ProCESS:::paths_to_sample(
      forest_data,
      highlight_sample
    )
    
    forest_data$highlight <-
      forest_data$to %in% highlight$to
  }
  
  
  
  edges <- forest_data %>%
    dplyr::select(
      "from",
      "to",
      "highlight"
    )
  
  
  
  graph <- tidygraph::as_tbl_graph(
    edges,
    directed = TRUE
  )
  
  
  
  graph <- graph %>%
    
    tidygraph::activate("nodes") %>%
    
    dplyr::left_join(
      
      forest_data %>%
        
        dplyr::rename(
          name = .data$to
        ) %>%
        
        dplyr::mutate(
          name = as.character(.data$name)
        ),
      
      by = "name"
    )
  
  
  
  graph <- graph %>%
    
    tidygraph::activate("edges") %>%
    
    dplyr::mutate(
      edge_alpha =
        tidygraph::.N()$alpha_level[to]
    )
  
  
  
  layout <- ggraph::create_layout(
    graph,
    layout = "tree",
    root = "WT"
  )
  
  
  
  max_Y <- max(
    layout$birth_time,
    na.rm = TRUE
  )
  
  layout$reversed_btime <-
    max_Y - layout$birth_time
  
  layout$y <-
    layout$reversed_btime
  
  
  
  samples_info <-
    forest$get_samples_info()
  
  nsamples <-
    nrow(samples_info)
  
  
  
  point_size <- c(
    0.5,
    rep(1, nsamples)
  )
  
  names(point_size) <- c(
    "N/A",
    samples_info %>%
      dplyr::pull(.data$name)
  )
  
  
  
  labels_every <-
    max_Y / 10
  
  
  
  tree_xlim <-
    range(layout$x, na.rm = TRUE) +
    c(-0.5, 0.5)
  
  
  
  
  graph_plot <- ggraph::ggraph(
    layout,
    "tree"
  ) +
    
    ggraph::geom_edge_link(
      
      edge_width = 0.1,
      
      ggplot2::aes(
        
        edge_color = ifelse(
          highlight,
          "black", #"indianred3",
          "black"
        ),
        
        alpha = .data$edge_alpha
      )
    )
  
  
  
  if (is.null(shape_label_function)) {
    
    graph_plot <- graph_plot +
      
      ggraph::geom_node_point(
        
        ggplot2::aes(
          color = .data$color_label,
          alpha = .data$alpha_level,
          size = .data$sample
        )
      )
    
  } else {
    
    graph_plot <- graph_plot +
      
      ggraph::geom_node_point(
        
        ggplot2::aes(
          color = .data$color_label,
          shape = .data$shape_label,
          alpha = .data$alpha_level,
          size = .data$sample
        )
      ) +
      
      ggplot2::scale_shape_manual(
        values = c(
          0:nsamples + 1
        )
      )
  }
  
  
  graph_plot <- graph_plot +
    
    ggplot2::scale_alpha_identity() +
    
    ggplot2::scale_color_manual(
      values = color_map
    ) +
    
    ggplot2::scale_size_manual(
      values = point_size
    ) +
    
    ggplot2::scale_x_continuous(
      limits = tree_xlim,
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      
      labels = seq(
        0,
        max_Y,
        labels_every
      ) %>%
        round() %>%
        rev(),
      
      breaks = seq(
        0,
        max_Y,
        labels_every
      ) %>%
        round()
    ) +
    
    ggplot2::theme_minimal() +
    
    ggplot2::theme(
      
      legend.position = "none",
      
      axis.line.x =
        ggplot2::element_blank(),
      
      axis.text.x =
        ggplot2::element_blank(),
      
      axis.ticks.x =
        ggplot2::element_blank(),
      
      plot.margin =
        ggplot2::margin(
          5, 5, 0, 5
        )
    ) +
    
    ggplot2::labs(
      color = NULL,
      shape = NULL,
      x = NULL,
      y = "Time"
    ) +
    
    ggplot2::guides(
      size = "none",
      shape = NULL,
      color = NULL
    )
  
  
  annotation_data <- layout %>%
    
    as.data.frame() %>%
    
    dplyr::filter(
      !is.na(.data$sample),
      .data$sample != "N/A"
    ) %>%
    
    dplyr::arrange(.data$x)
  
  
  annotation_plot_epi <- ggplot2::ggplot(
    
    annotation_data,
    
    ggplot2::aes(
      x = .data$x,
      y = 1,
      fill = .data$epistate
    )
    
  ) +
    
    ggplot2::geom_tile(
      width = 1,
      height = 1
    ) +
    
    ggplot2::scale_fill_manual(
      values = epistate_colors,
      drop = FALSE
    ) +
    
    ggplot2::scale_x_continuous(
      limits = tree_xlim,
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      expand = c(0, 0)
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        0, 5, 1, 5
      )
    )

  annotation_plot_geno <- ggplot2::ggplot(
    
    annotation_data,
    
    ggplot2::aes(
      x = .data$x,
      y = 1,
      fill = .data$mutant
    )
    
  ) +
    
    ggplot2::geom_tile(
      width = 1,
      height = 1
    ) +
    
    ggplot2::scale_fill_manual(
      values = mutant_colors,
      drop = FALSE
    ) +
    
    ggplot2::scale_x_continuous(
      limits = tree_xlim,
      expand = c(0, 0)
    ) +
    
    ggplot2::scale_y_continuous(
      expand = c(0, 0)
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(
        1, 5, 5, 5
      )
    )

  final_plot <-
    
    graph_plot /
    
    annotation_plot_epi /
    
    annotation_plot_geno +
    
    patchwork::plot_layout(
      heights = c(
        10,
        annotation_height,
        annotation_height
      )
    )
  
  
  return(final_plot)
}

plot_mutant_transition_graph <- function(
    sim,
    mutant_name,
    colors = phenotype_colors,
    edge_scale = 5
) {
  rates_df <- sim$get_rates()
  # Keep only switching events
  edges <- rates_df %>%
    filter(
      mutant == mutant_name,
      event == "switch"
    ) %>%
    transmute(
      from = epistate,
      to = first.child.epistate,
      weight = rate
    )
  
  # Make sure all phenotypes are present
  vertices <- data.frame(
    name = names(colors)
  )
  
  # Directed graph
  g <- graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = vertices
  )
  
  # Plot
  p <- ggraph(
    g,
    layout = "circle"
  ) +
    
    # fan separates reciprocal transitions
    geom_edge_fan(
      aes(
        width = weight,
        label = sprintf("%.2f", weight)
      ),
      arrow = arrow(
        length = unit(3, "mm"),
        type = "closed"
      ),
      start_cap = circle(7, "mm"),
      end_cap = circle(7, "mm"),
      colour = "grey40",
      label_colour = "black",
      label_size = 3.5,
      show.legend = FALSE
    ) +
    
    geom_node_point(
      aes(colour = name),
      size = 12
    ) +
    
    geom_node_text(
      aes(label = name),
      colour = "white",
      size = 4.5
    ) +
    
    scale_colour_manual(
      values = colors
    ) +
    
    # Same scale across G1/G2/G3
    scale_edge_width(
      limits = c(0, max(rates_df$rate[rates_df$event == "switch"])),
      range = c(0.5, edge_scale)
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = 0.35)
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = 0.35)
    ) +
    
    coord_fixed(
      clip = "off"
    ) +
    
    labs(
      title = mutant_name
    ) +
    
    theme_void() +
    
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}


df_to_phylo <- function(sample_forest, state_col = "epistate") {
  df <- sample_forest$get_nodes()
  # tips are cells that are never anyone's ancestor
  is_tip    <- !(df$cell_id %in% df$ancestor)
  tips      <- df$cell_id[is_tip]
  root      <- df$cell_id[is.na(df$ancestor)]
  internals <- c(root, setdiff(df$cell_id[!is_tip], root))  # root must be node n+1
  
  n   <- length(tips)
  ids <- setNames(c(seq_len(n), n + seq_along(internals)),
                  as.character(c(tips, internals)))
  
  ed <- df[!is.na(df$ancestor), ]
  bt <- setNames(df$birth_time, df$cell_id)
  
  phy <- list(
    edge        = unname(cbind(ids[as.character(ed$ancestor)],
                               ids[as.character(ed$cell_id)])),
    edge.length = unname(bt[as.character(ed$cell_id)] - bt[as.character(ed$ancestor)]),
    tip.label   = as.character(tips),
    node.label  = as.character(internals),
    Nnode       = length(internals)
  )
  class(phy) <- "phylo"
  
  phy <- reorder(phy)
  phy <- collapse.singles(phy)   # removes nodes that have only one child, adding up their branch lengths
  
  # tip states, in tip.label order
  s <- df[[state_col]][match(phy$tip.label, df$cell_id)]
  phy$states <- as.integer(factor(s))
  attr(phy$states, "levels") <- levels(factor(s))  # keeps the mapping, e.g. 1 = E1, 2 = E2
  phy
}
