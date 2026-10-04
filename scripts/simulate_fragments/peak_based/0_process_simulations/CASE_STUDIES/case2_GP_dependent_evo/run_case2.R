rm(list = ls())

library(ProCESS)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggalluvial)
library(patchwork)
library(scales)
library(grid)
library(ComplexHeatmap)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/CASE_STUDIES/case2_GP_dependent_evo")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/utils.R")
source("utils_case2.R")

set.seed(0)
dir.create("results", showWarnings = FALSE)
dir.create(file.path("results", "figures"), showWarnings = FALSE, recursive = TRUE)




# -----------------------------------------------------------------------------
# 1. Model definition
# -----------------------------------------------------------------------------
sim <- TissueSimulation(
  epigenetic_states = c("P1", "P2", "P3"),
  width = 1000, height = 1000,
  save_directory = TRUE,
  name = "case2"
)
sim$death_activation_level <- 50

rates_G1 <- list(
  P1 = list(duplication=.8, death=.1, P2=.8, P3=.7),
  P2 = list(duplication=.8, death=.1, P1=.7, P3=.8),
  P3 = list(duplication=.8, death=.1, P1=.8, P2=.7)
)
rates_G2 <- list(
  P1 = list(duplication=2, death=.1, P2=.02, P3=.02),
  P2 = list(duplication=2, death=.1, P1=.5, P3=.05),
  P3 = list(duplication=2, death=.1, P1=.5, P2=.05)
)
rates_G3 <- list(
  P1 = list(duplication=2.5, death=.1, P2=.05, P3=.5),
  P2 = list(duplication=2.5, death=.1, P1=.05, P3=.5),
  P3 = list(duplication=2.5, death=.1, P1=.02, P2=.02)
)

q_G1 <- get_epigenetic_Q(rates_G1, epistate_name="G1")
q_G2 <- get_epigenetic_Q(rates_G2, epistate_name="G2")
q_G3 <- get_epigenetic_Q(rates_G3, epistate_name="G3")

steady_states <- dplyr::bind_rows(
  as.data.frame(q_G1$steady_states) %>% mutate(clone="G1"),
  as.data.frame(q_G2$steady_states) %>% mutate(clone="G2"),
  as.data.frame(q_G3$steady_states) %>% mutate(clone="G3")
)
write.csv(steady_states, "results/steady_states.csv", row.names=FALSE)

# -----------------------------------------------------------------------------
# 2. Sequential clonal evolution
# -----------------------------------------------------------------------------
sim$add_mutant("G1", rates_G1)
sim$place_cell("G1[P2]", 500, 500)
sim$run_up_to_time(10)
p_t10 <- plot_tissue(sim) + facet_wrap(~epistate) + ggtitle("G1 population at t = 10")

sim$add_mutant("G2", rates_G2)
start_G2 <- sim$choose_border_cell_in("G1[P1]")
p_G2_origin <- plot_tissue(sim) +
  geom_point(data=start_G2, aes(x=position_x, y=position_y), inherit.aes=FALSE) +
  facet_wrap(~epistate) + ggtitle("Origin of G2")
sim$mutate_progeny(start_G2, "G2")
sim$set_rates(list("G1" = list(duplication=.4, death=.2)))
sim$run_up_to_time(40)
p_t40 <- plot_tissue(sim) + facet_wrap(~epistate) + ggtitle("After G2 expansion, t = 40")

sim$add_mutant("G3", rates_G3)
start_G3 <- sim$choose_border_cell_in("G2[P3]")
p_G3_origin <- plot_tissue(sim) +
  geom_point(data=start_G3, aes(x=position_x, y=position_y), inherit.aes=FALSE) +
  facet_wrap(~epistate) + ggtitle("Origin of G3")
sim$mutate_progeny(start_G3, "G3")
sim$set_rates(list("G2" = list(duplication=.4, death=.2)))
sim$set_rates(list("G1" = list(duplication=.2, death=.4)))
sim$run_up_to_time(60)
p_final_tissue <- plot_tissue(sim) + ggtitle("Final tissue, t = 60")

# -----------------------------------------------------------------------------
# 3. Spatial sample and forest
# -----------------------------------------------------------------------------
colors <- make_case2_color_maps()
bbox_width <- 40
bbox_p <- c(340, 290)
bbox_q <- bbox_p + bbox_width
p_sample_box <- plot_sample_box(sim, bbox_p, bbox_width)

sim$sample_cells("S3", bbox_p, bbox_q)
p_sample <- plot_tissue(sim, at_sample="S3") + ggtitle("Sample S3")
sample_forest <- sim$get_sample_forest()

sample_summary <- summarise_sample(sample_forest)
write.csv(sample_summary$epi_by_mutant, "results/sample_epistate_by_clone.csv", row.names=FALSE)
write.csv(sample_summary$mutant_by_epi, "results/sample_clone_by_epistate.csv", row.names=FALSE)

sankey_epi <- plot_epistate_alluvial(sample_summary$epi_by_mutant, colors$epi)
sankey_gen <- plot_mutant_alluvial(sample_summary$mutant_by_epi, colors$gen)

p_forest_epi <- plot_forest(sample_forest, color_map=colors$epistate) +
  theme(legend.position="none") + ggtitle("Sample forest: phenotype")
p_forest_gen <- plot_forest(sample_forest, color_map=colors$mutant) +
  theme(legend.position="none") + ggtitle("Sample forest: genetic clone")


p_evo <- wrap_plots(list=list(p_forest_gen,p_forest_epi,sankey_gen,sankey_epi),design = "aabb\naabb\nccdd") +
  plot_annotation(tag_levels = "a")


p_tran_mat <- wrap_plots(list=list(q_G1$plot_q,q_G2$plot_q,q_G3$plot_q))+
  plot_annotation(tag_levels = "a")
ggsave("results/figures/01_gen_phen_evolution.pdf", p_evo, width=12, height=8, dpi=300)
ggsave("results/figures/02_transition_mat.pdf", p_tran_mat, width=12, height=4, dpi=300)


sample_forest$save("results/sample_forest.sff")

# Objects used by the R Markdown report
save(sim, q_G1, q_G2, q_G3, steady_states, sample_forest, sample_summary,
     sankey_epi, sankey_gen, p_forest_epi, p_forest_gen, report,
     p_t10, p_G2_origin, p_t40, p_G3_origin, p_final_tissue,
     p_sample_box, p_sample,
     file="results/report_objects.RData")

# -----------------------------------------------------------------------------
# 4. Optional mutation placement
#    Keep credentials OUTSIDE the script. Set COSMIC_EMAIL and COSMIC_PASSWORD
#    in the environment before running this block.
# -----------------------------------------------------------------------------
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_references_v1.3.5/")

m_engine <- MutationEngine(
  setup_code="GRCh38", tumour_type="COADREAD", context_sampling=20,
  germline_subject="NA20514"
  # COSMIC_account=list(email=cosmic_email, password=cosmic_password)
)
mu_SNV <- 1e-9; mu_CNA <- 0; mu_INDELs <- 1e-9
CNA_G1_1 <- ProCESS::CNA(type="D", "11", from=48500001, len=1e7)
CNA_G1_2 <- ProCESS::CNA(type="D", "5", from=1000001, len=109499999, allele=0)
CNA_G1_3 <- ProCESS::CNA(type="D", "17", from=500001, len=42999999)
CNA_G2 <- ProCESS::CNA(type="D", "X", from=107000001, len=43499999)
CNA_G3_1 <- ProCESS::CNA(type="A", "8", from=127118340, len=1e7, src_allele=1)

passenger <- list(
  P1=c(SNV=mu_SNV, indel=mu_INDELs, CNA=mu_CNA),
  P2=c(SNV=mu_SNV, indel=mu_INDELs, CNA=mu_CNA),
  P3=c(SNV=mu_SNV, indel=mu_INDELs, CNA=mu_CNA)
)
m_engine$add_mutant("G1", passenger_rates=passenger,
                    drivers=list("APC Q1294Gfs*6", "TP53 R175H", CNA_G1_1, CNA_G1_2, CNA_G1_3))
m_engine$add_mutant("G2", passenger_rates=passenger,
                    drivers=list("KRAS Q22K", CNA_G2))
m_engine$add_mutant("G3", passenger_rates=passenger,
                    drivers=list("GNAS R844H", CNA_G3_1))
m_engine$add_exposure(time=0, coefficients=c(
  SBS1=.15, SBS5=.40, SBS18=.15, SBS17b=.20,
  ID1=.40, ID2=.40, ID18=.20, SBS88=.10
))
phylo_forest <- m_engine$place_mutations(
  sample_forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200
)
prj_dir <- "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/scATAC_project/ProCESS-scATAC/scripts/simulate_fragments/peak_based/0_process_simulations/process_GP_dependent_evo/"
setwd(prj_dir)
phylo_forest$save("results/phylo_forest.sff")

