### Analysis for DEgrade Targets - project month 1 ###
### Project lead: Christina Bligaard Pedersen
### April 2022

# Loading libraries
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(GGally)
library(reshape2)
library(cowplot)

# Set working directory
setwd('~/Documents/Intomics/')

# Read data set from Human Protein Atlas (RNA consensus tissue gene data)
expr_data <- read.table('rna_tissue_consensus.tsv', sep = '\t', header = T)

# Making a gene translation vector for translation between Ensembl ID and gene symbol becomes necessary
gene_trans <- expr_data %>%
  select(all_of(c('Gene', 'Gene.name'))) %>% 
  unique() %>% 
  tibble::deframe()

# Please note that gene symbols do not map one-to-one with the Ensembl IDs, which are all unique
# length(unique(gene_trans)) == length(gene_trans)
# length(unique(names(gene_trans))) == length(gene_trans)

# Reformatting data to a wide data frame (genes in rows, tissues in columns)
expr_df <- tidyr::spread(expr_data, key = Tissue, value = nTPM)

# Getting a list of tissues
tissues <- unique(expr_data$Tissue)

# Make the matrix
expr_mat <- expr_df %>%
  select(all_of(tissues)) %>%
  as.matrix()
rownames(expr_mat) <- expr_df$Gene



# Getting an idea about the data distribution to understand what constitutes high/low expression
ggplot(expr_data, aes(x = Tissue, y = nTPM)) +
  geom_boxplot() + 
  scale_y_continuous(trans = "log1p", breaks = c(0,100,1000,2000,5000,10000,100000)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + theme_bw()


### First, naive analysis (also good as sanity check) ----
# Simply calculating Pearson correlation to all genes
pcc <- sapply(names(gene_trans), function(gene) {cor(expr_mat[names(which(gene_trans == 'IGF2R')),], expr_mat[gene,], method = 'pearson')})

top_genes <- pcc[pcc > 0.75] %>% 
  sort(decreasing = T)

# length(top_genes) # 92 genes (incl. IGF2R itself)

# Translate to symbols and look at the top
gene_trans[names(top_genes)] %>% head()


# Plot an example scatter plots for top six non-self
plots <- lapply(2:7, function(i) {
  expr_mat %>% t() %>% as.data.frame() %>%
    ggplot(aes_string(x = names(top_genes[1]), y = names(top_genes[i]))) +
    geom_point() + theme_bw() +
    xlab(gene_trans[names(top_genes[1])]) + ylab(gene_trans[names(top_genes[i])]) +
    # ggtitle('Correlation across tissues') +
    scale_x_continuous(trans = 'log1p') + scale_y_continuous(trans = 'log1p')
})

plot_grid(plotlist = plots, ncol = 3)



