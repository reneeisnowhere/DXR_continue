#### a list of gene names was obtained using the
# search term "ZNF" from the website (HGNCname)https://www.genenames.org/
##on 1-20-2026


library(tidyverse)
library(readr)
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)



# data loading ------------------------------------------------------------

ZNF_search <- read_delim("C:/Users/renee/Documents/Ward Lab/DXR_project/ZNFhgnc-search-1768929420672_2026_01_20.txt",
                                                      delim = "\t", escape_double = FALSE,
                                                      trim_ws = TRUE)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")




# filtering out non-approved ----------------------------------------------



ZNF_list <- ZNF_search %>%
    dplyr::filter(status=="Approved") %>%
    dplyr::select(Symbol,ID)

gene_locs <- getBM(
  attributes = c("hgnc_symbol",
                 "ensembl_gene_id",
                 "chromosome_name",
                 "start_position",
                 "end_position",
                 "strand",
                 "gene_biotype"),
  filters = "hgnc_symbol",
  values  = ZNF_list$Symbol,
  mart    = ensembl)


# editing for non-standard and multiple rows ------------------------------


gene_locs <- gene_locs %>%
  filter(chromosome_name %in% c(1:22, "X", "Y")) %>%
  mutate(chr = paste0("chr", chromosome_name))
doubles <- gene_locs %>% count(hgnc_symbol) %>% dplyr::filter(n>1)
undouble <- gene_locs %>% dplyr::filter(hgnc_symbol %in% doubles$hgnc_symbol) %>%
  dplyr::filter(gene_biotype=="lncRNA")

filt_locs <- gene_locs %>%
  dplyr::filter(!ensembl_gene_id %in% undouble$ensembl_gene_id)
# convert to granges ------------------------------------------------------

gene_gr <- GRanges(seqnames = filt_locs$chr,
                   ranges=IRanges(start= filt_locs$start_position,
                                  end = filt_locs$end_position),
                   strand = ifelse(filt_locs$strand== 1, "+","-"),
                   hgnc_symbol=filt_locs$hgnc_symbol,
                   ensembl_gene_id=filt_locs$ensembl_gene_id,
                   gene_biotype = filt_locs$gene_biotype)

mcols(gene_gr)$name <- mcols(gene_gr)$hgnc_symbol
# export as bedfile -------------------------------------------------------


export(
  gene_gr,
  con = "C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/ZFP_proteins.bed",
  format = "BED")

