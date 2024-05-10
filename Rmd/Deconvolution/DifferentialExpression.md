Differential expression analysis
================
Marco Tello
2024-02-14

# Standard DEA

First we perform differential expression analysis utilizing only the
differences between diets to identify DEGs.

``` r
# Read expression data
exp_df <- fread(file = "~/GitHub/CellDecon/input/rnaseq/filtered_counts_Rnorv6.tsv")
rna_mat <- as.matrix(exp_df[,3:ncol(exp_df)])
rownames(rna_mat) <- exp_df$ENSEMBL

# Build metadata
condition <- str_replace(colnames(rna_mat), pattern = "_\\d+", replacement = "") 
coldat <- data.frame(Diet = as.factor(condition),
                     Sample = colnames(rna_mat), 
                     stringsAsFactors = TRUE)
coldat$Diet <- relevel(coldat$Diet, ref = "CSAA")
```

Use DESeq2 to identify differences between diets

``` r
deseq_rna_base <- DESeqDataSetFromMatrix(countData = rna_mat,
                                         colData = coldat,
                                         design = ~ Diet)

deseq_rna_base <- DESeq(deseq_rna_base)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
RNA_DE_PTS_vs_CSAA <- results(deseq_rna_base, 
                              name="Diet_PTS_vs_CSAA",
                              pAdjustMethod = "BH", 
                              alpha = 0.05)

RNA_DE_PTS_vs_CSAA <- cbind(rownames(RNA_DE_PTS_vs_CSAA),
                            RNA_DE_PTS_vs_CSAA)
names(RNA_DE_PTS_vs_CSAA)[1] <- "ENSEMBLID"
RNA_DE_PTS_vs_CSAA <- as.data.table(na.omit(RNA_DE_PTS_vs_CSAA))
RNA_DE_PTS_vs_CSAA[, padj := round(padj, 2)]
```

``` r
## Pull rat annotation in biomaRt
ensembl_v104 <- useEnsembl(biomart = "genes",
                           dataset = "rnorvegicus_gene_ensembl",
                           version = 104)

## Consider only the genes from the expression data
values <- unique(RNA_DE_PTS_vs_CSAA$ENSEMBLID)

##### One to one orthologs #####
conversion_table <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"),
                        filters = "ensembl_gene_id", 
                        values = values, 
                        mart= ensembl_v104)
setnames(conversion_table, "ensembl_gene_id", "ENSEMBL")
setnames(conversion_table, "external_gene_name", "GeneSymbol")

DEGs <- merge.data.table(x = conversion_table,
                         y =  RNA_DE_PTS_vs_CSAA, 
                         by.x = "ENSEMBL", 
                         by.y = "ENSEMBLID")
DEGs <- data.table(DEGs)

out_file <- "~/GitHub/CellDecon/output/DESeq2/base_DESeq_full.tsv"
if(!file.exists(out_file)){
  fwrite(DEGs, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}

out_file <- "~/GitHub/CellDecon/output/DESeq2/standard_DEG.tsv"
if(!file.exists(out_file)){
  fwrite(DEGs[padj <= 0.05], file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}
```

# Cell deconvolution results

First, we read the calculated cell fractions.

``` r
fractions <- fread("~/GitHub/CellDecon/output/CIBERSORTx/fractions_Sone2one_Mone2one_newExp.txt")
names(fractions) <- str_replace(str_replace(names(fractions), pattern = " Lineage", replacement = ""),pattern = " ", replacement = "")
fractions <- fractions[, .SD, .SDcols = c("Mixture","TCell","Macrophage","BCell","Neutrophil")]
fractions[, Mixture := gsub(pattern = "-", replacement = ".", x = Mixture)]
fractions <- data.table(fractions, key = "Mixture")

# Verify that samples are in the same order as the input rna_mat
fractions_sub <- fractions[colnames(rna_mat)]
fractions_norm <- fractions_sub[, .SD/sum(.SD), .SDcols = !"Mixture", by = Mixture]

# Format sample metadata 
coldat_decon <- merge.data.table(x = coldat, y = fractions_sub, 
                                 by.x = "Sample", by.y = "Mixture")
coldat_decon$Diet <- relevel(coldat_decon$Diet, ref = "CSAA")
coldat_decon <- data.table(coldat_decon, key = "Sample")
coldat_decon <- coldat_decon[colnames(rna_mat)]
```

Then we identify the differentially expressed genes incorporating the
cell fractions.

``` r
# Build DESeq object to perform the analysis
deseq_rna_decon <- DESeqDataSetFromMatrix(countData = rna_mat,
                                          colData = coldat_decon,
                                          design = ~ Diet + TCell + Macrophage + BCell)
deseq_rna_decon <- DESeq(deseq_rna_decon)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
decon_RNA_DE_PTS_vs_CSAA <- results(deseq_rna_decon, 
                                    name="Diet_PTS_vs_CSAA",
                                    pAdjustMethod = "BH", 
                                    alpha = 0.05)
decon_RNA_DE_PTS_vs_CSAA <- cbind(rownames(decon_RNA_DE_PTS_vs_CSAA),
                                  decon_RNA_DE_PTS_vs_CSAA)
names(decon_RNA_DE_PTS_vs_CSAA)[1] <- "Gene"
decon_RNA_DE_PTS_vs_CSAA <- as.data.table(na.omit(decon_RNA_DE_PTS_vs_CSAA))
```

``` r
## Pull rat annotation in biomaRt
ensembl_v104 <- useEnsembl(biomart = "genes",
                           dataset = "rnorvegicus_gene_ensembl",
                           version = 104)

## Consider only the genes from the expression data
values <- unique(decon_RNA_DE_PTS_vs_CSAA$Gene)

##### One to one orthologs #####
conversion_table <- getBM(attributes=c("ensembl_gene_id",
                                       "external_gene_name"),
                        filters = "ensembl_gene_id", 
                        values = values, 
                        mart= ensembl_v104)
setnames(conversion_table, "ensembl_gene_id", "ENSEMBL")
setnames(conversion_table, "external_gene_name", "Symbol")
conversion_table <- data.table(conversion_table)

decon_DEGs <- merge.data.table(x = conversion_table,
                               y =  decon_RNA_DE_PTS_vs_CSAA, 
                               by.x = "ENSEMBL", 
                               by.y = "Gene")

out_file <- "~/GitHub/CellDecon/output/DESeq2/decon_DESeq_full.tsv"
if(!file.exists(out_file)){
  fwrite(decon_DEGs, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}

decon_DEGs <- data.table(decon_DEGs)

out_file <- "~/GitHub/CellDecon/output/DESeq2/decon_DEG.tsv"
if(!file.exists(out_file)){
  fwrite(decon_DEGs[padj <= 0.05], file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}
```
