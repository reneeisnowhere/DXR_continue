### Code to pull promoter regions 1 kb +/- around TSS of KZFPs for compute matrix

library(GenomicRanges)
library(tidyverse)
library(rtracklayer)
library(biomaRt)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)



# ensembl <- useEnsembl(
#   biomart = "genes",
#   dataset = "hsapiens_gene_ensembl",
#   mirror = "www" )

KZFP_2026_list <- read_excel("data/Other_paper_data/human_protein_coding_KZFPs_Davis_natgen_2026.xlsx")
KZFP_2017_list <- read_excel("data/Other_paper_data/41586_2017_BFnature21683_MOESM103_ESM.xlsx")

human_KZFP <- KZFP_2026_list %>%
  dplyr::select(Label, `Gene ID`, chrom, start, end, strand) %>%
  dplyr::rename("Gene_ID"=`Gene ID`)
##biomart was down so I moved to pulling from Ens.Db
# results <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
#                                 "chromosome_name","transcript_start",
#                                 "transcript_end",
#                                 "strand",
#                                 "transcript_is_canonical"),
#                  filters = "ensembl_gene_id",
#                  values= human_KZFP$Gene_ID,
#                  mart= ensembl)
### retreive
edb <- EnsDb.Hsapiens.v86

tx <- transcripts(edb,filter = GeneIdFilter(human_KZFP$Gene_ID), return.type = "GRanges")
tx_by_gene <- split(tx,tx$gene_id)

tx_best <-do.call(c, lapply(tx_by_gene, function(x) {
  canon <- mcols(x)$tx_is_canonical

  canon <- if (is.null(canon)) {
    rep(0, length(x))
  } else {
    as.integer(ifelse(is.na(canon), 0, canon))
  }

  # rank manually
  score <- canon * 1e9 + width(x)

  x[which.max(score)]
}))

tx_best <- GRangesList(tx_best)
tx_best_gr <- unlist(tx_best, use.names = FALSE)
tss <- resize(tx_best_gr, width = 1, fix = "start")
export(tss, "data/Bed_exports/KZFP_tss_ensembl.bed")

