library(data.table)
library(biomaRt)
library(edgeR)
library(org.Rn.eg.db)


input_file <- "~/GitHub/CellDecon/input/rnaseq/filtered_counts_Rnorv6.tsv"
cibersort_path <- "~/GitHub/CellDecon/input/cibersortx"

##### Read shared data
# Read expression data
exp_df <- fread(input_file)
# Select only expression columns
exp_mat <- exp_df[, .SD, .SDcols = !c("ENSEMBL", 
                                      "GeneSymbol")]
# Perform CPM normalization on expression data
exp_mat <- cpm(exp_mat, log = FALSE)
# Merge in single DF
exp_df <- cbind(exp_df[, .SD, .SDcols =  c("GeneSymbol","ENSEMBL")], exp_mat)


# Read metaddata reference
phenotype_data <- fread(file = file.path(cibersort_path, 
                                         "Haemopedia-Mouse-RNASeq_samples.txt"))
# Select initial cell types to deconvolute
initial_celltypes <- sort(c("Neutrophil Lineage",
                            "T Cell Lineage",
                            "B Cell Lineage", 
                            "Macrophage Lineage"))
phenotype_sub <- phenotype_data[cell_lineage %in% initial_celltypes]

# Select only samples associated with the target lienages
reference_data <- fread(file = file.path(cibersort_path, 
                                         "Haemopedia-Mouse-RNASeq_cpm.txt"), 
                        select = c("GeneID", phenotype_sub$sampleId))

## Pull rat annotation in biomaRt
ensembl_v104 <- useEnsembl(biomart = "genes",
                           dataset = "rnorvegicus_gene_ensembl",
                           version = 104)

## Consider only the genes from the expression data
values <- unique(exp_df$GeneSymbol)

##### One to one orthologs #####
ortholog_1to1  <- getBM(attributes=c("ensembl_gene_id",
                                     "external_gene_name",
                                     "mmusculus_homolog_ensembl_gene",
                                     "mmusculus_homolog_associated_gene_name",
                                     "mmusculus_homolog_orthology_type",
                                     "mmusculus_homolog_orthology_confidence"),
                        filters = "external_gene_name", 
                        values = values, 
                        mart= ensembl_v104)
ortholog_1to1 <- as.data.table(ortholog_1to1)
ortholog_1to1 <- ortholog_1to1[mmusculus_homolog_orthology_type == "ortholog_one2one"
                                ][mmusculus_homolog_orthology_confidence == 1
                                  ][, .SD, 
                                    .SDcols = c("ensembl_gene_id",
                                                "external_gene_name",
                                                "mmusculus_homolog_ensembl_gene")]
out_file <- file.path(cibersort_path, "one2one_homologs.tsv")
if(!file.exists(out_file)){
  fwrite(ortholog_1to1, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}


### Subset data to orthologs #####
## Subset reference blood data
reference_sub <- merge.data.table(x = ortholog_1to1,
                                  y = reference_data, 
                                  by.x = "mmusculus_homolog_ensembl_gene", 
                                  by.y = "GeneID", 
                                  all = FALSE) 

out_file <- file.path(cibersort_path, "ortholog_reference_exp_cpm_one2one.tsv")
if(!file.exists(out_file)){
  fwrite(reference_sub, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}

cibersort_reference <- reference_sub[, mmusculus_homolog_ensembl_gene := NULL]
cibersort_reference <- cibersort_reference[, external_gene_name := NULL]
setnames(cibersort_reference, old = "ensembl_gene_id", new = "genes")
cibersort_reference <- na.omit(cibersort_reference)

out_file <- file.path(cibersort_path, "cibersortx_ortholog_ref_cpm_4cells_one2one.txt")
if(!file.exists(out_file)){
  fwrite(cibersort_reference, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}


## Subset expression data
# exp_df_sub <- merge.data.table(x = ortholog_1to1,
#                                y = exp_df, 
#                                by.x = "ensembl_gene_id", 
#                                by.y = "ENSEMBLID", 
#                                all = FALSE) 

cibersort_exp <- exp_df[, .SD, .SDcols = !c("GeneSymbol")] 
setnames(cibersort_exp, old = "ENSEMBL", "Gene")

out_file <- file.path(cibersort_path, "cibersortx_ortholog_mix_cpm_one2one.txt")
if(!file.exists(out_file)){
  fwrite(cibersort_exp, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}



##### Prepare extra files #####
# Create matrix with default value
pheno_matrix<-matrix(data = 2, 
                     nrow = length(unique(phenotype_sub$cell_lineage)),
                     ncol = length(phenotype_sub$sampleId))
# Assign celltypes to rownames
rownames(pheno_matrix) <- unique(phenotype_sub$cell_lineage)
sample_order <- names(cibersort_reference[, 
                                          .SD, 
                                          .SDcols = -c("genes")])

# Now we sort the phenotype data using the order of the reference data
phenotype_sub <- data.table(phenotype_sub, key = "sampleId")
phenotype_sub <- phenotype_sub[sample_order]

# Get what celltype correspond to each sample (simple parsing)
samples <- phenotype_sub$cell_lineage
# Subsitute default value in matrix in the corresponding column
for(cell_type in rownames(pheno_matrix)){
  col_index <- which(samples == cell_type)
  row_index <- which(rownames(pheno_matrix) == cell_type)
  pheno_matrix[row_index,col_index] <- 1
}
pheno_matrix <- data.table(initial_celltypes, pheno_matrix)

# Write pheno matrix
out_file <- file.path(cibersort_path, "cibersortx_ortholog_pheno_4cells.txt")
if(!file.exists(out_file)){
  fwrite(pheno_matrix, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = FALSE)
}

# ##### One to many orthologs #####
# ## Get ortholog attributes from the rat annotation only for genes in exp_df
# ortholog_table <- getBM(attributes=c("ensembl_gene_id",
#                                      "external_gene_name",
#                                      "mmusculus_homolog_ensembl_gene",
#                                      "mmusculus_homolog_associated_gene_name"),
#                         filters = "ensembl_gene_id", 
#                         values = values, 
#                         mart= ensembl_v110)
# ortholog_table <- as.data.table(ortholog_table)
# 
# ortholog_table <- ortholog_table[mmusculus_homolog_ensembl_gene != ""
#                                  ][external_gene_name != ""]
# 
# out_file <- file.path(cibersort_path, "one2many_homologs.tsv")
# if(!file.exists(out_file)){
#   fwrite(ortholog_table, file = out_file, 
#          append = FALSE, quote = FALSE, 
#          sep = '\t', na = "",
#          row.names = FALSE, col.names = TRUE)
# }
# 
# ### Subset data to orthologs #####
# ## Subset reference blood data
# reference_sub <- merge.data.table(x = ortholog_table[, 
#                                                      .SD, 
#                                                      .SDcols = c("ensembl_gene_id", 
#                                                                  "mmusculus_homolog_ensembl_gene")],
#                                   y = reference_data, 
#                                   by.x = "mmusculus_homolog_ensembl_gene", 
#                                   by.y = "GeneID", 
#                                   all = FALSE) 
# 
# out_file <- file.path(cibersort_path, "ortholog_reference_exp_cpm_one2many.tsv")
# if(!file.exists(out_file)){
#   fwrite(reference_sub, file = out_file, 
#          append = FALSE, quote = FALSE, 
#          sep = '\t', na = "",
#          row.names = FALSE, col.names = TRUE)
# }
# 
# cibersort_reference <- reference_sub[!duplicated(mmusculus_homolog_ensembl_gene)]
# cibersort_reference <- cibersort_reference[!duplicated(ensembl_gene_id)]
# cibersort_reference[, mmusculus_homolog_ensembl_gene := NULL]
# setnames(cibersort_reference, old = "ensembl_gene_id", new = "genes")
# cibersort_reference <- na.omit(cibersort_reference)
# 
# out_file <- file.path(cibersort_path, "cibersortx_ortholog_ref_cpm_one2many.txt")
# if(!file.exists(out_file)){
#   fwrite(cibersort_reference, file = out_file, 
#          append = FALSE, quote = FALSE, 
#          sep = '\t', na = "",
#          row.names = FALSE, col.names = TRUE)
# }
# 
# ## Subset expression data
# exp_df_sub <- merge.data.table(x = ortholog_table[, 
#                                                   .SD, 
#                                                   .SDcols = c("ensembl_gene_id", 
#                                                               "mmusculus_homolog_ensembl_gene")],
#                                y = exp_df, 
#                                by.x = "ensembl_gene_id", 
#                                by.y = "ENSEMBLID", 
#                                all = FALSE) 
# 
# cibersort_exp <- exp_df_sub[, .SD, .SDcols = !c("mmusculus_homolog_ensembl_gene")] 
# cibersort_exp <- cibersort_exp[!duplicated(ensembl_gene_id)]
# setnames(cibersort_exp, old = "ensembl_gene_id", "Gene")
# 
# out_file <- file.path(cibersort_path, "cibersortx_ortholog_mix_cpm_one2many.txt")
# if(!file.exists(out_file)){
#   fwrite(cibersort_exp, file = out_file, 
#          append = FALSE, quote = FALSE, 
#          sep = '\t', na = "",
#          row.names = FALSE, col.names = TRUE)
# }
