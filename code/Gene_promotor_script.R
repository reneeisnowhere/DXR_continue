### file to generate promoter regions around the human genome:
### i am only interested in generating one region around a gene, not the individual transcripts




library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(dplyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_gr <- genes(txdb)

## keeping standard chromosomes for genes

gene_gr <- keepSeqlevels(
  gene_gr,
  paste0("chr", c(1:22, "X", "Y")),
  pruning.mode = "coarse"
)
## getting a single promotor per transcript:

promotor_gr <- promoters(gene_gr, upstream = 2000, downstream = 500)

### keeping standard chromosomes

promotor_gr <- keepStandardChromosomes(promotor_gr, pruning.mode="coarse")

### keeping the standards chr1-22 and X and Y
keep_chr <- paste0("chr", c(1:22,"X","Y"))

promotor_gr <- keepSeqlevels(promotor_gr,intersect(seqlevels(promotor_gr),keep_chr),pruning.mode="coarse")

#### for specific promoter region around gene
export(promotor_gr,
       "data/hg38_gene_promotors_2kb_up_500bp_down.bed",format="BED")


gene_gr <- keepSeqlevels(gene_gr, intersect(seqlevels(gene_gr), keep_chr), pruning.mode="coarse")

#### for full gene body:
export(gene_gr,
       "data/hg38_genes_full.bed",format="BED")

# Make a 1-bp interval at the strand-aware TSS
tss_gr <- promoters(gene_gr,  upstream = 0,  downstream = 1)


# BED name and score fields
mcols(tss_gr)$name <- names(tss_gr)
mcols(tss_gr)$score <- 0


export(tss_gr,"hg38_TSS_1bp.bed", format = "BED")
head(tss_gr)
width(tss_gr) |> table()



