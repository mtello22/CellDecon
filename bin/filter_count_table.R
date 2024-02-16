library(data.table)

data_path <- "~/GitHub/CellDecon/input"
samples_CSAAPTS <- c("CSAA.1_024547.1", "CSAA.5_024548.5", 
                     "CSAA.9_024549.9", "CSAA.13_024550.13", 
                     "CSAA.17_024551.17", "CSAAPTS.3_024557.3",
                     "CSAAPTS.7_024558.7", "CSAAPTS.11_024559.11", 
                     "CSAAPTS.15_024560.15", "CSAAPTS.23_024561.23") 


# Load matrix with the feature counts
count_df <- fread(file.path(data_path, "rnaseq", "count_matrix_Rnorv6.tsv"))


# Remove columns with gene information 
exp_mat <- count_df[, .SD, .SDcols = !c("ENSEMBLID", "Length")]

# Select genes to keep
## Genes with at least min_reads in more than half of the samples 
min_reads <- 2
min_samples <- floor(ncol(exp_mat)/2)
genes_to_keep <- apply(exp_mat > min_reads, 
                       MARGIN = 1, 
                       sum) > min_samples

# Final gene count matrix
exp_mat <- exp_mat[genes_to_keep,]

exp_mat <- exp_mat[, .SD, .SDcols = samples_CSAAPTS]
exp_mat <- cbind(count_df[genes_to_keep, 1:2], exp_mat)


out_file <- file.path(data_path, "rnaseq", "filtered_counts.tsv")
if(!file.exists(out_file)){
  fwrite(exp_mat, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}

