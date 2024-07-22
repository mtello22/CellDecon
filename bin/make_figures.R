library(data.table)
library(ggplot2)

custom_volcano <- function(DEG_results, alpha, log2FC, pval_max = 16){
  num_na <- sum(is.na(DEG_results$padj))
  temp <- data.table(na.omit(DEG_results))
  temp[, alpha := ifelse(round(padj, 2) <= alpha, TRUE, FALSE)]
  temp[, log2FC := ifelse(abs(log2FoldChange) >= log2FC, TRUE, FALSE)]
  temp[, DEG := "No"]
  temp[, DEG := ifelse(alpha & !log2FC, "FDR", DEG)]
  temp[, DEG := ifelse(alpha & log2FC, "FDR and FC", DEG)]
  ggplot(temp, aes(x = log2FoldChange, y = -log10(padj), color = DEG)) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "gray", linewidth = 1) +
    geom_vline(xintercept = log2FC, linetype = "dashed", color = "gray", linewidth = 1) +
    geom_vline(xintercept = -log2FC, linetype = "dashed", color = "gray", linewidth = 1) +
    geom_point(alpha = 0.4, size = 3) +
    ylim(0, pval_max)+
    scale_x_continuous(breaks = seq(ceiling(range(temp$log2FoldChange))[1], ceiling(range(temp$log2FoldChange))[2])) +
    scale_color_manual(values = c("No" = "darkgray", "FDR" = "blue", "FDR and FC" = "red"), 
                       labels = c(paste("FDR (", as.character(table(temp$DEG)[1]), ")", sep = ""), 
                                  paste("FDR and FC (", as.character(table(temp$DEG)[2]),")", sep = ""), 
                                  paste("Not Sig. (", as.character(table(temp$DEG)[3] + num_na), ")", sep = ""))) +
    ylab(expression(-log[10]("adjusted p-value")))+
    xlab(expression(log[2]("fold change")))+
    labs(color = "DEG status") +
    theme_bw()
}


decon_deg <- fread("~/GitHub/CellDecon/output/DESeq2/decon_DESeq_full.tsv")
custom_volcano(decon_deg, 
               alpha = 0.05, 
               log2FC = 1)

ggplot(na.omit(decon_deg), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.05, breaks = seq(0, 1, by = 0.05), 
                 fill = "gray", color = "black") + 
  xlab(expression(italic("p-value")))+
  labs(title = "Histogram of p-values", y = "Frequency")+
  theme_bw()



base_degs <- fread("~/GitHub/CellDecon/output/DESeq2/base_DESeq_full.tsv")
custom_volcano(base_degs, 
               alpha = 0.05, 
               log2FC = 1, 
               pval_max = 2)


ggplot(na.omit(base_degs), aes(x = pvalue)) + 
  geom_histogram(binwidth = 0.05, breaks = seq(0, 1, by = 0.05), 
                 fill = "gray", color = "black") + 
  xlab(expression(italic("p-value")))+
  labs(title = "Histogram of p-values", y = "Frequency")+
  theme_bw()
