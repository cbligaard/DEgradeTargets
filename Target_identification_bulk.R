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
library(WGCNA)
library(hpar)


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
  scale_y_continuous(trans = "log1p", breaks = c(0,10,100,1000,10000,100000)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Note that most genes have expression in the range 0-100. The median is around 5.



### First, naive analysis (also good as sanity check) ----
# Simply calculating Pearson correlation from IGF2R to all genes (ignoring any NA's)
pcc <- sapply(names(gene_trans), function(gene) {stats::cor(expr_mat[names(which(gene_trans == 'IGF2R')),], expr_mat[gene,], method = 'pearson', use = "pairwise.complete.obs")})

top_genes <- pcc[pcc > 0.75] %>% 
  sort(decreasing = T)

# length(top_genes) # 92 genes (incl. IGF2R itself)

# Translate to symbols and look at the top
gene_trans[names(top_genes)] %>% head()


# Plot example scatter plots for top six non-self genes
plots <- lapply(2:7, function(i) {
  expr_mat %>% t() %>% as.data.frame() %>%
    ggplot(aes_string(x = names(top_genes[1]), y = names(top_genes[i]))) +
    geom_point() + theme_bw() +
    xlab(gene_trans[names(top_genes[1])]) + ylab(gene_trans[names(top_genes[i])]) +
    # ggtitle('Correlation across tissues') +
    scale_x_continuous(trans = 'log1p') + scale_y_continuous(trans = 'log1p')
})

plot_grid(plotlist = plots, ncol = 3)

# It looks good overall, also note that some of these other genes have generally higher expression than others





### A more well-thought out approach (based on WGCNA) ----

# All pairwise correlations have to be determined in order to generate networks
# Modules can then be defined by hierarchical clustering

# The following setting is needed according to WGCNA documentation
options(stringsAsFactors = FALSE)

# Check data quality
gsg <- goodSamplesGenes(t(expr_mat), verbose = 3)
table(gsg$goodGenes) # 225 genes are not "good"
table(gsg$goodSamples)

# We remove those genes and also transpose to fit WGCNA's functions
expr_mat_filtered <- t(expr_mat)[gsg$goodSamples, gsg$goodGenes]

# Check that IGF2R is still there
# names(which(gene_trans == 'IGF2R')) %in% colnames(expr_mat_filtered)


# To construct the network, we need to define a soft tresholding power suitable for 
# biological (typically scale free) networks

# Call the network topology analysis function with a set of powers
sft = pickSoftThreshold(expr_mat_filtered, powerVector = c(c(1:10), seq(from = 12, to=20, by=2)), verbose = 5)

# Plot the results
par(mfrow = c(1,2));

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")


# We pick the power 4 based on these plots and guidelines from WGCNA



# Now, we construct the network - largely using default values
# Analysis will be done in one block
net <- blockwiseModules(expr_mat_filtered, maxBlockSize = 21000,
                        power = 4, TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "DEgradeTx",
                        verbose = 3)

# Getting an overview of the total number of modules and their sizes
table(net$colors)


# Visualizing modules with dendrogram
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Saving some of the key information
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# Check which module IGF2R belongs to - and for sanity which modules, the top 20 correlated genes belong to
moduleLabels[names(which(gene_trans == 'IGF2R'))]
table(moduleLabels[names(top_genes[2:21])]) # Module 3 looks like the interesting one (1,882 genes)


# Identification of important genes in module 3
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expr_mat_filtered, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(expr_mat_filtered)));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

# Order on importance among genes in module
module3_membership <- geneModuleMembership[moduleLabels=="3",]
mod_order <- order(module3_membership$MM3, decreasing = T)

# Extract module 3 genes in order by importance
module3_genes <- rownames(module3_membership)[mod_order]

# Look for position of IGF2R
which(module3_genes == names(which(gene_trans == 'IGF2R'))) # Around third-in


# Make a table to show which genes co-express with IGF2R - simply taking all of module 3's genes here
cbind.data.frame("Gene" = module3_genes, "Gene.symbol" = gene_trans[module3_genes]) %>%
  write.csv(file = '~/Documents/Intomics/DEgradeTargets/Large_unfiltered_geneset.csv', row.names = F)



### Adding additional filtering steps on top of WGCNA ----

# Account for the median(?) expression level per gene to only get those with reasonably high expression (filter limit = 5(?))
median_expr_mod3 <- matrixStats::rowMedians(expr_mat[module3_genes,])

# Account for Pearson correlation to IGF2R (filter limit = 0.5(?))
pcc_mod3 <- pcc[module3_genes]


# Filter (somewhat crudely) and remake data tables
module3_subset1 <- cbind.data.frame("Gene" = module3_genes, 
                                    "Gene.symbol" = gene_trans[module3_genes], 
                                    "Median.expr" = median_expr_mod3,
                                    "PCC" = pcc_mod3) %>%
  filter(PCC > 0.5 & Median.expr > 5)

write.csv(module3_subset1, file = '~/Documents/Intomics/DEgradeTargets/Large_filtered_lvl1_geneset.csv', row.names = F)




# Get cellular location from HPA - we are only interested in extracellular (sectreted or transmembrane)
# Note: These data are highly complex, but it can support some filtering?

# Get data from subcellular location set
loc_mod3 <- getHpa(id = c(module3_genes),
                   hpadata = "hpaSubcellularLoc")
loc_relevant <- loc_mod3$Gene[grep('secreted|plasma membrane', tolower(loc_mod3$Main.location))]

# Get secretome data
sec_mod3 <- getHpa(id = c(module3_genes),
                   hpadata = "hpaSecretome")
sec_relevant <- sec_mod3$Gene[grep('secreted|membrane', tolower(sec_mod3$Protein.class))]


# Add location filtering and write table - very small set of genes kept here
module3_subset2 <- module3_subset1 %>%
  filter(Gene %in% loc_relevant & Gene %in% sec_relevant)

write.csv(module3_subset2, file = '~/Documents/Intomics/DEgradeTargets/Large_filtered_lvl2_geneset.csv', row.names = F)



# Make a plot for those passing the last filter
plots <- lapply(module3_subset2$Gene, function(g) {
  expr_mat %>% t() %>% as.data.frame() %>%
    ggplot(aes_string(x = names(top_genes[1]), y = g)) +
    geom_point() + theme_bw() +
    xlab(gene_trans[names(top_genes[1])]) + ylab(gene_trans[g]) +
    # ggtitle('Correlation across tissues') +
    scale_x_continuous(trans = 'log1p') + scale_y_continuous(trans = 'log1p')
})
plot_grid(plotlist = plots, ncol = 3)





### Look specifically for PCSK9, TARDBP, UCP2, DCN, APOD - plots and correlation ----
genes_of_interest <- c('PCSK9', 'TARDBP', 'UCP2', 'DCN', 'APOD')
plots <- lapply(genes_of_interest, function(g) {
  expr_mat %>% t() %>% as.data.frame() %>%
    ggplot(aes_string(x = names(top_genes[1]), y = names(which(gene_trans == g)))) +
    geom_point() + theme_bw() +
    xlab(gene_trans[names(top_genes[1])]) + ylab(g) +
    # ggtitle('Correlation across tissues') +
    scale_x_continuous(trans = 'log1p') + scale_y_continuous(trans = 'log1p')
})

plot_grid(plotlist = plots, ncol = 3)
pcc[names(gene_trans[match(genes_of_interest, gene_trans)])]




### Other data visualizations they can use ---
# Heatmaps?



### Future optimization: Do something with single cell data in high-expression tissue based on plots?
