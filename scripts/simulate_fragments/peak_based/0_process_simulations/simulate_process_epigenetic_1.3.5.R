rm(list=ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyverse)
# library(ComplexHeatmap)
# library(circlize)
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/plot_muller.R")
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC/utils/DLP.R")
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-DLP/scripts/utils/utils_plot.R")
# genetic_clades = read.csv(file = "RandomForest_PredictedClades.csv")
# genetic_clades_prop=genetic_clades %>% 
#   group_by(Clade_pred) %>% 
#   summarise(tot_cells=n()) %>% 
#   mutate(pct_cells=tot_cells/sum(tot_cells))
generate_shades <- function(color, n = 9) {
  stopifnot(n >= 2)

  # Create palette in Lab space for smooth color transitions
  pal <- grDevices::colorRampPalette(
    c("white",color),
    space = "Lab"
  )
  n_c=n+8
  x=pal(n_c)
  x=x[-c(1, length(x))]
  x=sample(x = x,size = n,replace = F)
  return(x)
}


# set the seed of the random number generator for repeatability
set.seed(0)
mutant_cols <- c(
  "0[E1]" = generate_shades("darkorange", 4)[1],
  "0[E2]" = generate_shades("darkorange", 4)[2],
  "0[E3]" = generate_shades("darkorange", 4)[3],
  "0[E4]" =generate_shades("darkorange", 4)[4],
  "A[E1]" = generate_shades("#ABB6FE", 4)[1],
  "A[E2]" = generate_shades("#ABB6FE", 4)[2],
  "A[E3]" = generate_shades("#ABB6FE", 4)[3],
  "A[E4]" =generate_shades("#ABB6FE", 4)[4],
  "B[E1]" = generate_shades("#465efdff", 4)[1],
  "B[E2]" = generate_shades("#465efdff", 4)[2],
  "B[E3]" = generate_shades("#465efdff", 4)[3],
  "B[E4]" =generate_shades("#465efdff", 4)[4],
  "C[E1]" = generate_shades("#0117a7ff", 4)[1],
  "C[E2]" = generate_shades("#0117a7ff", 4)[2],
  "C[E3]" = generate_shades("#0117a7ff", 4)[3],
  "C[E4]" =generate_shades("#0117a7ff", 4)[4],
  "D[E1]" = generate_shades("#9402eeff", 4)[1],
  "D[E2]" = generate_shades("#9402eeff", 4)[2],
  "D[E3]" = generate_shades("#9402eeff", 4)[3],
  "D[E4]" =generate_shades("#9402eeff", 4)[4]
)


mutant_cols <- c(
  "0[E1]" = "darkorange",
  "0[E2]" = "darkorange",
  "0[E3]" = "darkorange",
  "0[E4]" ="darkorange",
  "A[E1]" = "#ABB6FE",
  "A[E2]" = "#ABB6FE",
  "A[E3]" = "#ABB6FE",
  "A[E4]" ="#ABB6FE",
  "B[E1]" = "#465efdff",
  "B[E2]" = "#465efdff",
  "B[E3]" = "#465efdff",
  "B[E4]" ="#465efdff", 
  "C[E1]" = "#0117a7ff",
  "C[E2]" = "#0117a7ff",
  "C[E3]" = "#0117a7ff",
  "C[E4]" ="#0117a7ff",
  "D[E1]" = "#9402eeff",
  "D[E2]" = "#9402eeff",
  "D[E3]" = "#9402eeff",
  "D[E4]" ="#9402eeff"
)
# unlink("P05/")
# create a tissue simulation with two epigenetic states
sim <- TissueSimulation(epigenetic_states = c("E1", "E2","E3","E4"),
                        width = 1000, height = 1000,save_snapshots = F,name = "P05")

rates_genetic_clone_0 <- list(
  E1 = list(duplication = 0.6, death = 0.1),
  E2 = list(duplication = 0.6, death = 0.1),
  E3 = list(duplication = 0.6, death = 0.1),
  E4 = list(duplication = 0.6, death = 0.1)
)


# add a mutant "A" and set its species rates
rates_genetic_clone_A <- list(
  E1 = list(duplication = 0.8, death = 0.1, E2 = 0.3, E3 = 0.5, E4 = 0.01),
  E2 = list(duplication = 0.8, death = 0.1, E1 = 0.12, E3 = 0.3, E4 = 0.2),
  E3 = list(duplication = 0.8, death = 0.1, E1 = 0.01, E2 = 0.2, E4 = 0.02),
  E4 = list(duplication = 0.8, death = 0.1, E1 = 0.02, E2 = 0.2, E3 = 0.5)
)

rates_genetic_clone_B <- list(
  E1 = list(duplication = 1, death = 0.1, E2 = 0.3, E3 = 0.5, E4 = 0.1),
  E2 = list(duplication = 1, death = 0.1, E1 = 0.8, E3 = 0.5, E4 = 0.2),
  E3 = list(duplication = 1, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 1, death = 0.1, E1 = 0.2, E2 = 0.5, E3 = 0.02)
)

rates_genetic_clone_C <- list(
  E1 = list(duplication = 1.5, death = 0.1, E2 = 0.3, E3 = 0.2, E4 = 0.1),
  E2 = list(duplication = 1.5, death = 0.1, E1 = 0.8, E3 = 0.1, E4 = 0.1),
  E3 = list(duplication = 1.5, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 1.5, death = 0.1, E1 = 0.4, E2 = 0.5, E3 = 0.02)
)

rates_genetic_clone_D <- list(
  E1 = list(duplication = 2, death = 0.1, E2 = 0.3, E3 = 0.2, E4 = 0.1),
  E2 = list(duplication = 2, death = 0.1, E1 = 0.8, E3 = 0.1, E4 = 0.1),
  E3 = list(duplication = 2, death = 0.1, E1 = 0.4, E2 = 0.2, E4 = 0.1),
  E4 = list(duplication = 2, death = 0.1, E1 = 0.4, E2 = 0.5, E3 = 0.02)
)


sim$add_mutant("0", rates_genetic_clone_0)

# states <- names(rates_genetic_clone_B)
# 
# # Initialize matrix
# mat <- matrix(0,
#               nrow = length(states),
#               ncol = length(states),
#               dimnames = list(states, states))
# 
# # Fill in switching rates
# for (from in states) {
#   for (to in names(rates_genetic_clone_B[[from]])) {
#     if (!(to %in% c("duplication", "death"))) {
#       mat[from, to] <- rates_genetic_clone_B[[from]][[to]]
#     }
#   }
# }
# diag(mat) <- -rowSums(mat - diag(diag(mat)))
# 
# 
# 
# Heatmap(mat,col = colorRampPalette(c("white", "steelblue"))(10),cluster_rows = F,cluster_columns = F,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# place one cell of "A" in epigenetic state "E1"
sim$place_cell("0[E1]", 500, 500)

# let the simulation evolve until the species "A[E2]" has less than 10 cells
sim$run_up_to_time(20)


sim$add_mutant("A",rates_genetic_clone_A)
starting_cell=sim$choose_border_cell_in("0")
plot_tissue(sim)+
  # geom_vline(xintercept = starting_cell$position_x)+
  # geom_hline(yintercept = starting_cell$position_y)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)

sim$mutate_progeny(starting_cell,"A")
sim$set_rates(list("0" = list(duplication = 0.4)))
sim$run_up_to_time(40)
plot_tissue(sim)



# Sampling ncells with random box sampling of boxes of size n_w x n_h


sim$add_mutant("B",rates_genetic_clone_B)
starting_cell=sim$choose_border_cell_in("0")
plot_tissue(sim)+
  # geom_vline(xintercept = starting_cell$position_x)+
  # geom_hline(yintercept = starting_cell$position_y)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)

sim$mutate_progeny(starting_cell,"B")
sim$set_rates(list("A" = list(duplication = 0.4, death = 0.2)))
sim$set_rates(list("0" = list(duplication = 0.05, death = 0.8)))
sim$run_up_to_time(60)
plot_tissue(sim)


bbox_width=20
bbox1_p <- c(550, 500)
bbox1_q <- bbox1_p + bbox_width

library(ggplot2)

# view the boxes
plot_tissue(sim) +
  geom_rect(
    xmin = bbox1_p[1],
    xmax = bbox1_q[1],
    ymin = bbox1_p[2],
    ymax = bbox1_q[2],
    fill = NA,
    color = "black"
  )

sim$add_mutant("C",rates_genetic_clone_C)
starting_cell=sim$choose_border_cell_in("B",c(550, 500), c(570, 520))
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)

sim$mutate_progeny(starting_cell,"C")

sim$set_rates(list("A" = list(duplication = 0.2, death = 0.5)))
sim$set_rates(list("B" = list(duplication = 0.4, death = 0.5)))
sim$run_up_to_time(80)
plot_tissue(sim)+
  geom_point(data = starting_cell,aes(x=position_x,y=position_y),inherit.aes = FALSE)


sim$add_mutant("D",rates_genetic_clone_D)
sim$mutate_progeny(sim$choose_border_cell_in("C"),"D")
# sim$set_rates(list("A" = list(duplication = 0.2, death = .9)))
sim$set_rates(list("B" = list(duplication = 0.1, death = .5)))
sim$set_rates(list("C" = list(duplication = 0.1, death = .5)))
sim$run_up_to_time(100)
plot_tissue(sim,color_map = mutant_cols)
# plot_muller(sim,color_map = mutant_cols)


# genetic_clades_prop





# Sample Clone A
n_w <- n_h <- 15
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("A" = ncells), n_w, n_h)
sim$sample_cells("S_A", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim,at_sample = "S_A")

# Sample Clone B
n_w <- n_h <- 5
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("B" = ncells), n_w, n_h)
sim$sample_cells("S_B", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim,at_sample = "S_B")

# Sample Clone C
n_w <- n_h <- 40
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("C" = ncells), n_w, n_h)
sim$sample_cells("S_C", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim,at_sample = "S_C")

# Sample Clone D
n_w <- n_h <- 25
ncells <- 0.8 * n_w * n_h
bbox <- sim$search_sample(c("D" = ncells), n_w, n_h)
sim$sample_cells("S_D", bbox$lower_corner, bbox$upper_corner)
plot_tissue(sim,at_sample = "S_D")

sample_forest <- sim$get_sample_forest()




# sampled_genetic_prop= sample_forest$get_nodes() %>% 
#   filter(!is.na(sample)) %>% 
#   group_by(mutant) %>% 
#   summarise(n=n())  %>% 
#   mutate(pct_cells_sim=n/sum(n))
# 
# sampled_genetic_prop%>% 
#   ggplot(aes(x = "", y = n, fill = mutant)) +
#   geom_col(width = 1) +
#   scale_fill_manual(values = c("A"="#cf02eeff","B"="#e6ecf8ff","C"="#abb6feff","D"="#0117a7ff"))+
#   coord_polar(theta = "y") +
#   theme_void() +
#   labs(fill = "Mutant")
# 
# 
# genetic_clades_prop= genetic_clades_prop %>% 
#   mutate(mutant=case_when(Clade_pred==1~"A",
#                           Clade_pred==2~"B",
#                           Clade_pred==3~"C",
#                           TRUE~"D"))
# merged_df = inner_join(sampled_genetic_prop,genetic_clades_prop)
# merged_df %>% 
#   ggplot(aes(x=pct_cells,y = pct_cells_sim,colour = mutant))+
#   geom_point()+
#   scale_color_manual(values = c("A"="#cf02eeff","B"="#e6ecf8ff","C"="#abb6feff","D"="#0117a7ff"))+
#   theme_minimal()+
#   labs(x="Simulted data",y="Ground truth")
# 
# # plot it
plot_forest(sample_forest,color_map = mutant_cols) %>%
  annotate_forest(sample_forest)

# sample_forest_old=load_sample_forest("sample_forest_atac_epigenome_new_version_chloe_data.sff")
# plot_forest(sample_forest_old,color_map = mutant_cols) %>% 
#   annotate_forest(sample_forest_old)

sample_forest$save("sample_forest_atac_epigenome_1.3.5_pat05.sff")
dir.create("process_references_v1.3.5")
setwd("process_references_v1.3.5")

m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", context_sampling = 20,
                           germline_subject = "NA20514",
                           COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it","password"="2*db!XQ4sgQ!dbg"))


mu_SNV = 1e-9
mu_CNA = 0
mu_INDELs = 1e-9

CNA_Clone0_1 = ProCESS::CNA(type = "D", "11",
                            from = 48500001, len = 1e7)
CNA_Clone0_2  = ProCESS::CNA(type = "D", "5",
                             from = 1000001, len = 109499999,allele = 0)
CNA_Clone0_3 = ProCESS::CNA(type = "D","17",
                            from = 500001, len = 42999999)

# chrX      5 107000001 150500000       0.25         1 chrX:107000001:150500000 43499999 D

CNA_CloneA = ProCESS::CNA(type = "D", "X",
                          from = 107000001, len = 43499999)



CNA_CloneB_1 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_2 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_3 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_4 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_5 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_6 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_7 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)
CNA_CloneB_8 = ProCESS::CNA(type = "A", "8",
                            from = 127118340, len = 1e7,src_allele = 1)

CNA_CloneB_9 = ProCESS::CNA(type = "A", "8",
                            from = 26500001, len = 16499999)


CNA_CloneB_10 = ProCESS::CNA(type = "A", "7",
                            from = 130500001, len = 28499999)

CNA_CloneB_11 = ProCESS::CNA(type = "A", "15",
                             from = 88000001, len = 13999999)


CNA_CloneC = ProCESS::CNA(type = "D", "8",
                          from = 107000001, len = 22999999)

CNA_CloneD = ProCESS::CNA(type = "D", "14",
                          from = 46500001, len = 58999999)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "0",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(list("APC R1450*", allele = 1),CNA_Clone0_1,CNA_Clone0_2,CNA_Clone0_3))
m_engine$add_mutant(mutant_name = "A",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_CloneA))

m_engine$add_mutant(mutant_name = "B",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_CloneB_1,CNA_CloneB_2,CNA_CloneB_3,CNA_CloneB_4,CNA_CloneB_5,CNA_CloneB_6,CNA_CloneB_7,CNA_CloneB_8))

m_engine$add_mutant(mutant_name = "C",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_CloneC))


m_engine$add_mutant(mutant_name = "D",
                    passenger_rates = list("E1" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E2" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E3" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA),
                                           "E4" = c(SNV = mu_SNV, indel = mu_INDELs,CNA=mu_CNA)),
                    drivers = list(CNA_CloneD))
m_engine$add_exposure(time = 0,coefficients = c(SBS1 = 0.15,SBS5 = 0.40,
                                                SBS18 = 0.15,SBS17b = 0.20,ID1 = 0.40,ID2 = 0.40,ID18=0.2,SBS88 = 0.10))
phylo_forest <- m_engine$place_mutations(sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save("../phylo_forest_atac_epigenome_1.3.5_pat05.sff")



# # phylo_forest$get_bulk_allelic_fragmentation("S1") %>% View()
# phylo_forest$get_bulk_allelic_fragmentation() %>%
#   mutate(karyotype=paste(major,minor,sep=":")) %>% 
#   filter(karyotype!="1:1") %>% 
#   mutate(seg_id=paste(chr,begin,end,sep=":")) %>% 
#   mutate(tot_cn=major+minor) %>% 
#   ggplot(aes(y=seg_id,x=ratio,fill=as.factor(tot_cn)))+
#   geom_col()+
#   labs(x="CCF",y="Segment",fill="CN")+
#   theme_minimal()
# 
# 
# 
# single_cell_cnas = phylo_forest$get_cell_allelic_fragmentation() %>% 
#   mutate(karyotype=paste(major,minor,sep=":")) %>% 
#   filter(karyotype!="1:1")
# sample_forest$get_nodes() %>% 
#   group_by(mutant,epistate) %>% 
#   count() %>% 
#   group_by(mutant) %>% 
#   mutate(tot_cell=sum(n)) %>% 
#   mutate(pct_epi=n/tot_cell) %>% 
#   ggplot(aes(x=mutant,y=pct_epi,fill=epistate))+
#   geom_bar(stat="identity")
# 
# 
