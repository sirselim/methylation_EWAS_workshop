#### methylation_EWAS_workshop
## R script for exploring EWAS analysis of whole blood and buccal cells
## Author: Miles Benton (m.benton@qut.edu.au)
## Date created: 3rd May 2017
## Date modified: 4th May 2017
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## this workshop uses public Illumina 450K data from several studies 
## (GSE40279, GSE48472)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#####################################################
### 1. load required packages
#####################################################
# if you need to install uncomment the steps below
# source("https://bioconductor.org/biocLite.R")
# biocLite("minfi")
require(minfi)

#####################################################
### 2. load data
#####################################################
meth.anno <- readRDS('data/Ilmn450K_anno.RDS')
betas <- readRDS('data/blood_buccal.RDS')
# load('data/blood_buccal.RData')

#####################################################
### 3. set up samples/data and basic QC explore
#####################################################

# define samples
tissues <- c(rep(0, 5), rep(1, 5))
tissues <- as.factor(tissues)
levels(tissues) <- c('Blood', 'Buccal')

# create a densoty plot
densityPlot(betas, sampGroups = tissues)

# Hierarchical Clustering
d <- t(betas)
rownames(d) <- tissues
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on all probes", labels = tissues)

#####################################################
### 4. run dmpfinder - logistic regression
#####################################################

# minfi dmpfinder
results <- dmpFinder(betas, pheno = tissues, type = "categorical", qCutoff = 1)
results <- na.omit(results)
head(results)
tail(results)

# explore some results
table(results$qval < 0.00005)
top.results <- results[results$qval <= 0.00005,]
tail(top.results)

# bonferroni correction
table(results$pval <= (0.05/450000))
# top.results <- results[results$pval <= (0.05/450000),]

#####################################################
### 5. plotting the results
#####################################################

## very basic plotting 
# top CpG marker ('best'/most significantly associated)
plotCpg(betas, cpg = "cg05350268", pheno = tissues, type = "categorical", measure = "beta", ylab = "methylation (beta)", ylim = c(0,1))

# bottom CpG marker ('best'/most significantly associated)
plotCpg(betas, cpg = "cg07309136", pheno = tissues, type = "categorical", measure = "beta", ylab = "methylation (beta)", ylim = c(0,1))

## boxplots
# lets look at a few boxplots
boxplot(betas['cg05350268',] ~ tissues, col = c('darkred', 'cadetblue'), ylim = c(0,1))
boxplot(betas['cg19653345',] ~ tissues, col = c('darkred', 'cadetblue'), ylim = c(0,1))
boxplot(betas['cg03467405',] ~ tissues, col = c('darkred', 'cadetblue'), ylim = c(0,1))
boxplot(betas['cg07046426',] ~ tissues, col = c('darkred', 'cadetblue'), ylim = c(0,1))

## multiple boxplots
# loop through a few sites for the top markers (most significantly associated)
tissue.sites <- rownames(top.results)[c(1:9)]

par(mfrow = c(3,3))
for (cpg in tissue.sites) {
  
  boxplot(betas[cpg,] ~ tissues, col = c('darkred', 'cadetblue'), ylab = "methylation (beta)", xlab = "tissue", main = cpg, ylim = c(0,1))
  
}

## multiple boxplots
# loop through a few sites for the bottom markers (least significantly associated)
tissue.sites <- rownames(tail(results, n = 12))[c(1:9)]

par(mfrow = c(3,3))
for (cpg in tissue.sites) {
  
  boxplot(betas[cpg,] ~ tissues, col = c('darkred', 'cadetblue'), ylab = "methylation (beta)", xlab = "tissue", main = cpg, ylim = c(0,1))
  
}

## multiple boxplots
# can also write this out to a file (pdf)
pdf('figures/blood_v_buccal_topsites_plot.pdf', width = 7, height = 10)
tissue.sites <- rownames(top.results)[c(1:12)]

par(mfrow = c(4,3))
for (cpg in tissue.sites) {
  
  boxplot(betas[cpg,] ~ tissues, col = c('darkred', 'cadetblue'), ylab = "methylation (beta)", xlab = "tissue", main = cpg)
  
}
dev.off()

## Hierarchical Clustering 2
# subset the beta matrix by the top markers and cluster
sig.betas <- betas[rownames(betas) %in% rownames(top.results),]

# Hierarchical Clustering 
d <- t(sig.betas)
rownames(d) <- tissues
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on selected top tissue discrimination CpG sites", labels = tissues)

#####################################################
### 6. Prepare and run pathways analysis
#####################################################

## create a list of genes which underlie the CpG sites extracted
# explore annotation
head(meth.anno) # see what the annotation looks like

# subset annotation based on the top.results
anno.sub <- meth.anno[meth.anno$IlmnID %in% rownames(top.results),]

# how many different features can be explored?
table(anno.sub$CHR) # number of significant sites per chromosome
length(unique(unlist(strsplit(anno.sub$UCSC_RefGene_Name, split = ';')))) # number of unique genes

# get the gene list from the provided 450K annotation object
genes <- unique(unlist(strsplit(anno.sub$UCSC_RefGene_Name, split = ';')))

# write list of genes out to a file
write.table(genes, 'EWAS_blood_buccal_genelist.txt', row.names = F, col.names = F, quote = F)

# this list can now be uploaded to your pathways analysis tool
# suggest WEB-based GEne SeT AnaLysis Toolkit (http://webgestalt.org)

####/END