###This is documentation on how I got the KZFP regions defined for my signal calculations across each KZFP
## I pulled the Supplementary 2 table from:
##Imbeault, M., Helleboid, PY. & Trono, D. KRAB zinc-finger proteins contribute to the evolution of gene regulatory networks.
#Nature 543, 550–554 (2017). https://doi.org/10.1038/nature21683

library(tidyverse)
library(readxl)
library(GenomicRanges)
library(rtracklayer)

KZFP_2017_list <- read_excel("data/Other_paper_data/41586_2017_BFnature21683_MOESM103_ESM.xlsx")
View(KZFP_2017_list)


df_out <- KZFP_2017_list %>%
  dplyr::filter(Species == "Homo_sapiens.GRCh38") %>%
  dplyr::select(Label, chrom, start,end, strand) %>%
  mutate(chrom =paste0("chr",chrom))

gr_out <-   GRanges(
  seqnames = df_out$chrom,
  ranges = IRanges(start = df_out$start, end = df_out$end),
  strand = df_out$strand,
  name = df_out$Label)

export(gr_out, "data/Bed_exports/KZFP.bed", format = "BED")
