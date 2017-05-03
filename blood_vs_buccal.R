####
##

require(glmnet)
require(FactoMineR)
require(caret)
require(magrittr)
require(ROCR)


#
load('data/blood_buccal.RData')



#####################################################
### Bootstrapping glmnet
#####################################################


## setting up for training
x <- t(na.omit(betas))

# outcome for tissue compare
y <- c(rep(0, 5), rep(1, 5))
y <- as.factor(y)
levels(y) <- c('Blood', 'Buccal')


# Do glmnet
fit <- glmnet(x, y, family = "binomial", alpha = 0.9)
# Get model coefficients for glmnet
Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
# Get CpG list for which coefficients are not 0
cpgs <- rownames(betas[Coefficients@i,])

plot(fit)
plot(fit, xvar="lambda")

# create a results object for the above
result <- betas[Coefficients@i,]
nrow(result)  # number of snps extracted 

# define a snp set for PCA 
pca.data <- t(betas[Coefficients@i,])
pca.data <- as.data.frame(pca.data)

# add additional covars to plot in PCA 
pca.data$cell_type <- y

# define the covars for PCA 
quali.covars <- c(grep('cell_type', colnames(pca.data)))

# run PCA
res.pca <- PCA(pca.data, quali.sup = quali.covars, graph = FALSE)
# plot
pdf('figures/blood_bucal_glmnet_0.9_pca_161025.pdf', paper = 'a4')
# par(cex = 0.65)
# plot(res.pca, habillage = ncol(pca.data), col.hab = c("darkgreen","cadetblue", "darkred"))
# plotellipses(res.pca, axis = c(1,2), cex = 0.65, col.hab = c("darkgreen","cadetblue", "darkred"), label = 'ind.sup')
plotellipses(res.pca, keepvar = c("cell_type"), label = 'none')
# look at individual covars
# plotellipses(res.pca, keepvar = "IQ", label = 'ind.sup')
# plotellipses(res.pca, keepvar = "sex", label = 'ind.sup')
# plotellipses(res.pca, keepvar = "etnicity", label = 'ind.sup')
dev.off()
####


boxplot(betas[grep('cg24380059', rownames(betas)),] ~ y, main = 'cg24380059')


pdf('figures/blood_buccal_differences_160210.pdf', width = 7, height = 10)
cpg.list <- cpgs
par(mfrow = c(3,3), las = 2, cex.main = 0.85)
for (i in cpg.list) {
  
  boxplot(betas[grep(i, rownames(betas)),] ~ y, ylim = c(0, 1), col = c('cadetblue', 'darkred'), 
          main = paste(i, "\n", paste(unique(unlist(strsplit(meth.anno[grep(i, meth.anno$IlmnID),]$UCSC_RefGene_Name, split = ';'))), collapse = '; ')))
  grid()
  
}
dev.off()


unique(unlist(strsplit(meth.anno[grep(i, meth.anno$IlmnID),]$UCSC_RefGene_Name, split = ';')))

bb.markers <- read.csv('blood_vs_buccal_bootNet_1000iter_160210.csv', head = T)
bb.markers.top <- bb.markers[bb.markers$iterations == 1000,]$IlmnID

pdf('figures/blood_buccal_differences_bootNet1000_160210.pdf', width = 7, height = 10)
cpg.list <- bb.markers.top
par(mfrow = c(3,3), las = 2, cex.main = 0.85)
for (i in cpg.list) {
  
  boxplot(betas[grep(i, rownames(betas)),] ~ y, ylim = c(0, 1), col = c('cadetblue', 'darkred'), 
          main = paste(i, "\n", paste(unique(unlist(strsplit(meth.anno[grep(i, meth.anno$IlmnID),]$UCSC_RefGene_Name, split = ';'))), collapse = '; ')))
  grid()
  
}
dev.off()


###################################################################
###### on methead
###################################################################

# packages
require(marmalaid)
require(glmnet)
require(caret)
require(magrittr)
library(foreach)
library(doParallel)

# setting up parallel processing
# cl <- makeCluster(2)
# registerDoParallel(cl)
# foreach(i=1:3) %dopar% sqrt(i)
registerDoParallel(cores = 3)

# 450k anno
meth.anno <- read.csv('HumanMethylation450_15017482_v1-2.csv', head = T, as.is = T, skip = 7)

# marmalaid anno
data(annotation_v1.1)

# WB from population exp
GSE40279 <- annotation[grep('GSE40279', annotation$GSE),]
GSE40279.betas <- getbeta(samples=GSE40279$Id, probes=meth.anno$IlmnID, marmalaiddir='/publicdata/marmalaid/')

# buccal
buccal <- annotation[grep('Buccal', annotation$TISSUE),]
buccal.betas <- getbeta(samples=buccal$Id[c(1:25, 30:75)], probes=meth.anno$IlmnID, marmalaiddir='/publicdata/marmalaid/')

## setting up for training
beta.matrix <- na.omit(cbind(GSE40279.betas, buccal.betas))
x <- t(beta.matrix)

# outcome for tissue compare
y <- c(rep(0, 656), rep(1, 71))
y <- as.factor(y)
levels(y) <- c('Blood', 'Buccal')

## bootstrap function
bootNet <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample){
  # load packages
  require(glmnet)
  # require(survival)
  cpg_list <- list()
  # bootstrap process 
  for (i in 1:iter){
    set.seed(i)
    # Select a random number of patients
    # this is currently only set up for 30 samples (15 control vs 15 case), need to fix this
    # newDataInd <- c(sample(1:15,floor(sub_sample*15)), sample(16:30,floor(sub_sample*15))) # sampling from both groups (quantitative)
    newDataInd <- c(sample(grep('Blood', outcome), floor(sub_sample*(length(grep('Blood', outcome))/2))), 
                    sample(grep('Buccal', outcome), floor(sub_sample*(length(grep('Buccal', outcome))/2))))
    #Subset the data
    newData <- data[newDataInd,]
    # In the outcome variable get the same patients as were selected for this iteration
    newOut <- outcome[newDataInd]
    # Do glmnet
    fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    # Get model coefficients for glmnet
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- rownames(beta_matrix[Coefficients@i,])
    name <- paste('run:', i, sep = '')
    cpg_list[[name]] <- cpgs
    print(i)
  }
  return(cpg_list)
  cat('\n', ' ...Processing Done...')
}

# run bootstrapping (100 iterations)
#
cpgs.bootstrap <- bootNet(data = x, outcome = y, Alpha = 0.1, iter = 10, beta_matrix = beta.matrix, sub_sample = 0.666)

# generate a data frame that can then be annotated
bootData <- data.frame(IlmnID = names(sort(table(unlist(cpgs.bootstrap)))),  iterations = sort(table(unlist(cpgs.bootstrap))))
bootData <- bootData[order(-bootData$iterations),]
table(bootData$iterations == 10)
table(bootData$iterations >= 9)

# pre tissue
markers.bloodbuccal <- bootData
table(markers.bloodbuccal$iterations >= 80)
markers.bloodbuccal.top <- markers.bloodbuccal[markers.bloodbuccal$iteration == 10,]


pdf('Big_blood_buccal_differences_bootNet1000_160210.pdf', width = 7, height = 10)
cpg.list <- markers.bloodbuccal.top$IlmnID
par(mfrow = c(3,3), las = 2, cex.main = 0.85)
for (i in cpg.list) {
  
  boxplot(beta.matrix[grep(i, rownames(beta.matrix)),] ~ y, ylim = c(0, 1), col = c('cadetblue', 'darkred'), 
          main = paste(i, "\n", paste(unique(unlist(strsplit(meth.anno[grep(i, meth.anno$IlmnID),]$UCSC_RefGene_Name, split = ';'))), collapse = '; ')))
  grid()
  
}
dev.off()

########
## parallel bootstrap function
########
bootNet.par <- function(data, outcome, Alpha, iter, Lambda, beta_matrix, sub_sample, cores){
  # load packages
  require(glmnet)
  library(foreach)
  library(doParallel)
  # register cores
  registerDoParallel(cores = cores)
  # bootstrap process 
  foreach (i = 1:iter, .combine = c) %dopar% {
    set.seed(i)
    # Select a random number of patients/samples
    newDataInd <- c(sample(grep('Blood', outcome), floor(sub_sample*(length(grep('Blood', outcome))/2))), 
                    sample(grep('Buccal', outcome), floor(sub_sample*(length(grep('Buccal', outcome))/2))))
    #Subset the data
    newData <- data[newDataInd,]
    # In the outcome variable get the same patients as were selected for this iteration
    newOut <- outcome[newDataInd]
    # Do glmnet
    fit <- glmnet(x = newData, y = newOut, family = "binomial", alpha = Alpha)
    # Get model coefficients for glmnet
    Coefficients <- coef(fit, s = 0.001)  # if cv is performed this can be coef(fit, s = cv.fit$lambda.min)
    # Get CpG list for which coefficients are not 0
    cpgs <- rownames(beta_matrix[Coefficients@i,])
    list(cpgs)
  }
}

# run bootstrapping (100 iterations)
# methead has 12 cores, so try with 10 for now
system.time(
  cpgs.bootstrap <- bootNet.par(data = x, outcome = y, Alpha = 0.1, iter = 2500, beta_matrix = beta.matrix, sub_sample = 0.666, cores = 10)
)
# this took 5hrs on methead
# WARNING: uses a LOT of RAM!! Got up to ~100GB when run across 10 cores

# generate a data frame that can then be annotated
bootData <- data.frame(IlmnID = names(sort(table(unlist(cpgs.bootstrap)))),  iterations = sort(table(unlist(cpgs.bootstrap))))
bootData <- bootData[order(-bootData$iterations),]
table(bootData$iterations == 2500)
# table(bootData$iterations >= )

pdf('Big_blood_buccal_differences_bootNet2500_160216.pdf', width = 7, height = 10)
cpg.list <- bootData[bootData$iterations == 2500,]$IlmnID
par(mfrow = c(3,3), las = 2, cex.main = 0.85)
for (i in cpg.list) {
  
  boxplot(beta.matrix[grep(i, rownames(beta.matrix)),] ~ y, ylim = c(0, 1), col = c('cadetblue', 'darkred'), 
          main = paste(i, "\n", paste(unique(unlist(strsplit(meth.anno[grep(i, meth.anno$IlmnID),]$UCSC_RefGene_Name, split = ';'))), collapse = '; ')))
  grid()
  
}
dev.off()

###################################################################
###### /END/ on methead
###################################################################











