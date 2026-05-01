## Concatenation of all KZFP hits using fimo
library(data.table)

files <-  list.files(
  path = "C:/Users/renee/Other_projects_data/DXR_data/H3K27ac_KZFP_fimo",
  pattern = "best_site\\.narrowPeak$",
  recursive = TRUE,
  full.names = TRUE
)
dt <- rbindlist(lapply(files, fread))


fwrite(dt, "C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/KZFP_all_best_sites.narrowPeak", sep = "\t")


system("wc -l chunk_*/best_site.narrowPeak")
system("wc -l all_best_sites.narrowPeak")
length(files)
sapply(files, function(f) file.info(f)$size)

KZFP_motif_info <- readRDS("data/RDS_files/KZFP_motif_info.RDS")
df <- read.table("C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/KZFP_fixed.narrowPeak", header = FALSE, sep = "\t")
colnames(df)[4] <- "Motif_ID"

df_annot <- merge(df, KZFP_motif_info, by="Motif_ID", all.x=TRUE)

df_annot$motif_label <- paste0(
  df_annot$Motif_ID,
  ":",
  df_annot$TF_Name)

df_short <-
df_annot %>%
  dplyr::select("V1","V2","V3","motif_label","V5","V6","V7","V8")
split_list <- split(df_short, df_short$motif_label)

dir.create("C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/KZFP_by_motif", showWarnings = FALSE)
safe_names <- gsub("[:/\\*?\"<>|]", "_", names(split_list))
for (i in seq_along(split_list)) {
  write.table(
    split_list[[i]][, 1:8],
    file = file.path(
      "C:/Users/renee/R_STUFF/DXR_continue/data/Bed_exports/KZFP_by_motif",
      paste0(safe_names[i], ".narrowPeak")
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}


sapply(split_list, nrow)
head(split_list$`M00272_3.00:ZNF200`)
