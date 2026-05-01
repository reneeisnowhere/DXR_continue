### Conversion of bed to fasta file

library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)


bed <- import("C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/repeatmasker_sets/ERV1/MER61C_ERV1_LTR.bed", format = "BED")


genome <- BSgenome.Hsapiens.UCSC.hg38
seqs <- getSeq(genome, bed)

names(seqs) <- paste0(seqnames(bed), ":", start(bed), "-", end(bed))

names(seqs) <- mcols(bed)$name

writeXStringSet(seqs, filepath="C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/fasta_files/MER61C_ERV1_LTR.fa")
