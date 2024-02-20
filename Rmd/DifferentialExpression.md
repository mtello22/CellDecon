Differential expression analysis
================
Marco Tello
2024-02-14

``` r
conversion_table <- fread("~/GitHub/CellDecon/input/rnaseq/conversion_table.tsv")
```

# Standard DEA

First we perform differential expression analysis utilizing only the
differences between diets to identify DEGs.

``` r
# Read expression data
exp_df <- fread(file = "~/GitHub/CellDecon/input/rnaseq/filtered_counts.tsv")
rna_mat <- as.matrix(exp_df[,3:ncol(exp_df)])
rownames(rna_mat) <- exp_df$ENSEMBLID

# Build metadata
condition <- str_replace(colnames(rna_mat), pattern = ".\\d+", replacement = "") 
condition <- str_replace(condition, pattern = "_\\d{6}.\\d+", replacement = "") 
condition <- str_replace(condition, pattern = "CSAAPTS", replacement = "PTS") 

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
names(RNA_DE_PTS_vs_CSAA)[1] <- "Gene"
RNA_DE_PTS_vs_CSAA <- as.data.table(na.omit(RNA_DE_PTS_vs_CSAA))
RNA_DE_PTS_vs_CSAA[, padj := round(padj, 2)]

out_file <- "~/GitHub/CellDecon/output/DESeq2/base_DESeq_full.tsv"
if(!file.exists(out_file)){
  fwrite(RNA_DE_PTS_vs_CSAA, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}
```

## Visualize results

``` r
hist(RNA_DE_PTS_vs_CSAA$pvalue, main = "CSAA vs PTS", xlab = "p-value")
```

![](DifferentialExpression_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
alpha = 0.05
log2FC = 1
custom_volcano(RNA_DE_PTS_vs_CSAA, alpha, log2FC)
```

![](DifferentialExpression_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
base_DEGs <- merge.data.table(x = conversion_table,
                              y =  RNA_DE_PTS_vs_CSAA, 
                              by.x = "ENSEMBL", 
                              by.y = "Gene")

base_DEGs[, padj := round(padj, 2)]

out_file <- "~/GitHub/CellDecon/output/DESeq2/standard_DEG.tsv"
if(!file.exists(out_file)){
  fwrite(base_DEGs[padj <= 0.05], file = out_file, 
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
fractions <- fractions[, .SD, .SDcols = c("Mixture","BCell","TCell","Macrophage","Neutrophil")]
fractions[, Mixture := gsub(pattern = "-", replacement = ".", x = Mixture)]
fractions <- data.table(fractions, key = "Mixture")

# Verify that samples are in the same order as the input rna_mat
fractions_sub <- fractions[colnames(rna_mat)]
fractions_norm <- fractions_sub[, .SD/sum(.SD), .SDcols = !"Mixture", by = Mixture]

# Format sample metadata 
coldat_decon <- merge.data.table(x = coldat, y = fractions_sub, 
                                 by.x = "Sample", by.y = "Mixture")
coldat_decon$Diet <- relevel(coldat_decon$Diet, ref = "CSAA")
```

Then we identify the differentially expressed genes incorporating the
cell fractions.

``` r
# Build DESeq object to perform the analysis
deseq_rna_decon <- DESeqDataSetFromMatrix(countData = rna_mat,
                                          colData = coldat_decon,
                                          design = ~ Diet + TCell + Macrophage + BCell )
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

out_file <- "~/GitHub/CellDecon/output/DESeq2/decon_DESeq_full.tsv"
if(!file.exists(out_file)){
  fwrite(decon_RNA_DE_PTS_vs_CSAA, file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}
```

## Visualize deconvolution results

``` r
hist(decon_RNA_DE_PTS_vs_CSAA$pvalue, main = "CSAA vs PTS", xlab = "p-value")
```

![](DifferentialExpression_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
alpha = 0.05
log2FC = 1
custom_volcano(decon_RNA_DE_PTS_vs_CSAA, alpha, log2FC)
```

![](DifferentialExpression_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
decon_DEGs <- merge.data.table(conversion_table, 
                               decon_RNA_DE_PTS_vs_CSAA, 
                               by.x = "ENSEMBL", 
                               by.y = "Gene")

decon_DEGs[, padj := round(padj, 2)]

out_file <- "~/GitHub/CellDecon/output/DESeq2/decon_DEG.tsv"
if(!file.exists(out_file)){
  fwrite(decon_DEGs[padj < 0.05], file = out_file, 
         append = FALSE, quote = FALSE, 
         sep = '\t', na = "",
         row.names = FALSE, col.names = TRUE)
}
```