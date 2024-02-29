# Rat Blood Cell Deconvolution Project

This repository contains the code for the cell deconvolution project. The goal is to identify differentially expressed genes between rats fed with a regular diet vs PTS supplemented diet.

Data to reproduce the results is not available.

Under the **\~/Rmd** folder, I organized the reports necessary to reproduce the results:

-   **DataExploration/EDA.md**: Report for the exploratory data analysis of the RNA-seq data

-   **Deconvolution/DifferentialExpression.md**: Report that performs differential expression analysis with/without estimated cell fractions using DESeq2. To visualize DEGs, I included a volcano plot and a histogram of p-values.

-   **DataExploration/CIBERSORTxSignatureMatrix.md**: Report that visualizes the signature matrix and cell proportions estimated by CIBERSORTx. It also compares the weight of DEGs in the signature matrix to determine if important genes for deconvolution were also identified as DEGs.

-   **Interpretation/fGSEA.md**: Report that performs fast Gene Set Enrichment Analysis in the DESeq2 deconvoluted results against Reactome Pathways.
