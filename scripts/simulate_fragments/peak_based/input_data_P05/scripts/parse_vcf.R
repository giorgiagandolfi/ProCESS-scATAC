library(vcfR)

vcf=read.vcfR(file = 'input_data_P05/data/HIPEC_P05_N123_consensus_PASS_vep_filtered.vcf')
vcf@fix

tb = vcfR::vcfR2tidy(vcf,info_only = T)



fix_field = tb[["fix"]] %>%
  dplyr::rename(
    chr = CHROM,
    from = POS,
    ref = REF,
    alt = ALT) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    from = as.numeric(from),
    to = from + nchar(alt)) %>%
  dplyr::ungroup() %>%
  dplyr::select(chr, from, to, ref, alt, dplyr::everything())

vep_field = tb[['meta']] %>%
  dplyr::filter(ID == "CSQ") %>%
  dplyr::select(Description) %>%
  dplyr::pull()


tmp_vep_field = strsplit(vep_field, split = "|", fixed = TRUE) %>% unlist()
vep_field = tmp_vep_field[1:length(tmp_vep_field)-1]

# Tranform the fix field by splittig the CSQ and select the columns needed
fix_field= fix_field %>%
  dplyr::mutate(CSQ = strsplit(CSQ, ",")) %>%
  tidyr::unnest(CSQ) %>%
  tidyr::separate(CSQ, vep_field, sep = "\\|",convert = T) %>%
  dplyr::select(chr, from, to, ref, alt, IMPACT, SYMBOL, Gene, dplyr::everything(),-VariantKey)  #can add other thing, CSQ, HGSP

fix_field_position = fix_field %>% 
  select(chr,from,to,ref,alt,SYMBOL,HGVSc,IMPACT
         ) %>% 
  distinct()
 fix_field_position
