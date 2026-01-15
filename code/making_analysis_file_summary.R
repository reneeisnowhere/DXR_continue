### This is for combining all QC data for a summary of all samples in one file

library(tidyverse)

alignment_summary <- read_delim("data/alignment_summary.txt")
all_peak_final <- read_delim("data/all_peak_final_dataframe.txt",delim = "\t")
frag_peak_summary <- read_delim("data/number_frag_peaks_summary.txt")
sample_info <- read_delim("data/sample_info.tsv")
MQC_summary <- read_delim("data/multiqc_data_trim/Summary_of_multiqcfiles.txt")
peak_ct <- read_delim("data/peaks/peaks_cts_FINAL.txt", delim = "\t")

# test <-
sample_info %>%
  left_join(.,all_peak_final, by=join_by(Histone_Mark, Individual, Treatment, Timepoint) ) %>%
  left_join(., alignment_summary,by = join_by(Histone_Mark, Individual, Treatment, Timepoint)) %>%
  dplyr::filter(Treatment != "5FU") %>%
  left_join(., peak_ct, by=c("Library ID"="Sample")) %>%
  left_join(., frag_peak_summary, by = c("Library ID"="sample")) %>%
  write_delim(., "data/full_summary_QC_metrics.txt",delim = "\t")


MQC_summary %>%
  dplyr::select(prefix, Histone_Mark,Individual,Treatment, Timepoint,`FastQC_mqc-generalstats-fastqc-total_sequences`) %>%
  distinct()
