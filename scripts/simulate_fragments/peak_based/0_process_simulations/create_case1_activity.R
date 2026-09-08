library(dplyr)
library(tidyr)

set.seed(123)

# your full list of hallmark pathways (example — replace with your actual list)
activity <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway.rds")
activity_df <- convert_activity_list(activity_list = activity)
pathways <- activity_df$pathway %>% unique()

epistates <- c("E1", "E2", "E3")

# assign each epistate a distinct, non-overlapping set of "signature" pathways
n_signature <- 10
set.seed(123)
shuffled <- sample(pathways)

signature_list <- list(
  E1 = shuffled[1:n_signature],
  E2 = shuffled[(n_signature+1):(2*n_signature)],
  E3 = shuffled[(2*n_signature+1):(3*n_signature)]
)

# build the dataframe: high activity for signature pathways, low/baseline for the rest
df <- expand_grid(epistate = epistates, pathway = pathways) %>%
  rowwise() %>%
  mutate(
    is_signature = pathway %in% signature_list[[epistate]],
    activity = if (is_signature) {
      rnorm(1, mean = 0.85, sd = 0.05)   # high activity for that epistate's pathways
    } else {
      rnorm(1, mean = 0.30, sd = 0.08)   # low/baseline activity otherwise
    }
  ) %>%
  ungroup() %>%
  mutate(activity = pmin(pmax(activity, 0), 1)) %>%  # clip to [0,1]
  select(epistate, pathway, activity)
result <- split(df, df$epistate) %>%
  lapply(function(x) setNames(x$activity, x$pathway))
saveRDS(object = result,file = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/reference_data/a_scores_pathway_case1.rds")
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

wide <- df %>%
  pivot_wider(names_from = epistate, values_from = activity) %>%
  column_to_rownames("pathway")

mat <- as.matrix(wide)

col_fun <- colorRamp2(c(0, 1), c("white", "steelblue3"))
Heatmap(t(mat), col = col_fun, name = "activity",column_names_gp = gpar(fontsize = 8))
