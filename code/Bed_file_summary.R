library(rtracklayer)
library(dplyr)
library(purrr)

# Choose directory interactively
bed_dir <- rstudioapi::selectDirectory(caption = "Select folder with BED files")

# Get all BED files
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

# Function to summarize one BED file
summarize_bed <- function(file) {
  gr <- import(file, format = "BED")
  widths <- width(gr)

  tibble(
    filename = basename(file),
    n_lines = length(gr),
    average_length = mean(widths),
    max_len = max(widths),
    min_len = min(widths)
  )
}

# Apply across files
summary_df <- map_dfr(bed_files, summarize_bed)

# View + save
print(summary_df)
write.csv(summary_df, file.path("data/Other_paper_data/ERV1_file_summary.csv"), row.names = FALSE)
