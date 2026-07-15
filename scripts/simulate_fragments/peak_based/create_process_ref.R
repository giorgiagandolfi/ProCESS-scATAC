library(ProCESS)
library(dplyr)
set.seed(12345)
outdir <- 'process_references'
setwd(outdir)
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "COADREAD", 
                           context_sampling = 20,COSMIC_account = list("email"="giorgia.gandolfi@phd.units.it",
                                                                       "password"="2*db!XQ4sgQ!dbg"))
