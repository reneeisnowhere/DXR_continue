### How I converted from hg37 to hg38 coordinates
# Sample GSM2466657		Query DataSets for GSM2466657
# Status	Public on Mar 06, 2017
# Title	ZNF780A
# Sample type	SRA
#
# GSE78099	ChIP-exo of human KRAB-ZNFs transduced in HEK 293T cells and KAP1 in hES H1 cells
# Source name	transduced 293T cell line
# Organism	Homo sapiens
# Characteristics	host cell line: HEK 293T
# tagged protein: ZNF780A

library(rtracklayer)
library(GenomicRanges)
library(liftOver)
#read in file
bed <- import("C:/Users/renee/Other_projects_data/DXR_data/ZNF780A_peak_data/GSM2466657_ZNF780A_peaks_processed_score_signal_exo.bed.gz", format = "BED")
 ##download hg19to hg38.over.change

chain <- import.chain("C:/Users/renee/Other_projects_data/DXR_data/ZNF780A_peak_data/hg19ToHg38.over.chain")



lifted <- liftOver(bed, chain)

table(elementNROWS(lifted))
mapped_idx <- elementNROWS(lifted) > 0
bed_mapped <- bed[mapped_idx]
lifted_gr <- unlist(lifted[mapped_idx])

length(lifted_gr)
export(lifted_gr, "C:/Users/renee/R_STUFF/DXR_continue/data/Other_paper_data/GSM2466657_lifted_hg38_peaks_ZNF780A.bed")
