### Set up ###
rm(list = ls())
#setwd("C:/Users/Kalvin Edwin/Desktop/STA426/Project")


#### Load required libraries ####
library(MASS)
library(tidyverse)
library(edgeR)
library(limma)
library(devtools)
library(markerGeneProfile)
library(dtangle)
library(glmnet)
library(MineICA)
library(amritr)
library(mixOmics)
library(preprocessCore)
#library(declassification) # could not find this package
library(DSA)


### Read Data FGCZ ###
info <- read.table('PreFrontalCortex_infos.tsv', sep = '\t', header = TRUE)
nrow(info) # 102 patients

data <- read.table('Count_QC-raw-count.txt', sep = '\t', header = TRUE)
nrow(data) # rows: counts recorded from 21'505 genes, columns: count reads for each patient

# library(Rtsne)
# set.seed(10)
# data2 <- data[1:1000,]
# tsne <- Rtsne(data2, check_duplicates = FALSE, pca = FALSE, perplexity = 40, theta = 0.5)
# plot(tsne$Y[,1], tsne$Y[,2])


#### Data Preprocessing ####

# Create a new matrix with unique genes by removing the duplicated reads (here line index). 
# The first time I tried to make the gene_name as rownames, it did not work because some reads were duplicated.
# It automatically gave me the name of those genes, so I knew which ones to remove.
# I removed the read with the lower expression values
data_unique <- data[-c(8382, 5347, 21408, 21409, 21406, 3914, 21412, 21412, 21402, 21403, 21410, 20201, 21399,10000,21404,
                        21463, 5918, 21407, 3036, 8486, 21398, 7989, 21400, 2706, 21401, 21405, 6682, 7324, 21461, 1982,
                        21051, 21462, 21411, 3706),]

# Create a new matrix with gene_name as row names
rownames(data_unique) <- data_unique[,2]
data_unique <- data_unique[,3:ncol(data_unique)]

# Filter low expressed genes
group <- c(rep(1,52),rep(2,50))
dge <- DGEList(counts = data_unique, group = group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalize
geneExpr <- cpm(dge)
geneExpr <- as.data.frame(geneExpr)
geneResid = geneExpr - rowMeans(geneExpr)


#### Read first marker gene ####

# Load Darmanis dataset
load("GSE67835.RData")

# Keep common genes
int = intersect(rownames(geneExpr),rownames(Darmanis))
y <- DGEList(counts= Darmanis[int,])
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Darmanis[,1]/sum(Darmanis[,1])),0.75))
Darmanis <- cpm(y, log=FALSE)

DarmanisMean = sapply(unique(DarmanisCells),function(x)rowMeans(Darmanis[,names(which(DarmanisCells==x))]))

int = intersect(rownames(geneExpr),rownames(DarmanisMean))
design = model.matrix(~.-1,data.frame(cell= DarmanisCells[colnames(Darmanis)]))
colnames(design) = gsub('cell','',colnames(design))

v <- voom(Darmanis[int,rownames(design)], design, plot=FALSE)
fit <- lmFit(v, design)
x = sapply(colnames(design),function(x)paste(x,'-(',paste(colnames(design)[colnames(design)!=x],collapse = "+"),')/4',sep = ''))
con = makeContrasts(contrasts = x,levels = colnames(design))
fit = contrasts.fit(fit,con)
fit <- eBayes(fit, robust=TRUE)

DarmanisMarkers = NULL
for(i in 1:ncol(design)){
  x = fit$coef[p.adjust(fit$p.v[,i],'fdr')<0.05,i]
  DarmanisMarkers[[colnames(design)[i]]] =  names(head(sort(-x[x>0]),100))
}
unlist(lapply(DarmanisMarkers,length))
DarmanisMarkersFull = NULL
for(i in 1:ncol(design)){
  x = fit$coef[p.adjust(fit$p.v[,i],'fdr')<0.05,i]
  DarmanisMarkersFull[[colnames(design)[i]]] =  names(sort(-x[x>0]))
}
unlist(lapply(DarmanisMarkersFull,length))


#### Deconvolution ####
sce = DarmanisMean
ge = geneExpr

commongenes <- intersect (rownames(ge), rownames(sce))
ge <- ge[pmatch(commongenes, rownames(ge)), ]
sce <- sce[pmatch(commongenes, rownames(sce)), ]

y <- cbind(sce, ge) # Normally it's 2^ge but it leads to 'inf' values, so maybe something is missing in the preprocess part
#y[sapply(y, simplify = 'matrix', is.infinite)] <- 1e100 # should we apply a threshold for 'inf' value ?

y <- DGEList(counts=y)
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Darmanis[,1]/sum(Darmanis[,1])),0.75))
y <- cpm(y, log=FALSE)

#y = normalizeQuantiles(y)

pure_samples = as.list(1:5)
names(pure_samples) = colnames(y)[1:5]
markers = lapply(DarmanisMarkers,intersect,rownames(y))

# dtangle
dtDarmanis <- dtangle(log2(t(y)), pure_samples=pure_samples, markers = markers[names(pure_samples)])$estimates
#matplot(markers,dtDarmanis[6:nrow(dtDarmanis),], xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

# DSA
#dsaDarmanis <- getcompM(geneResid,method="DSA",marker=markers) # inbuilt function unknown
dsaDarmanis <- Deconvolution(y, )

# NNLS
dat = t(y[unlist(markers),])
lsDarmanis = apply(dat[-c(1:5),],1,function(x)coef(glmnet(t(dat[1:5,unlist(markers)]),x,lambda = 0, lower.limits = 0,intercept = FALSE,standardize = FALSE))[-1,])
lsDarmanis = t(lsDarmanis)/colSums(lsDarmanis)

# CIBERSORT
Z = y[unlist(DarmanisMarkers),-c(1:5)]
X = y[unlist(DarmanisMarkers),1:5]
X <- (X - mean(X))/sd(as.vector(X))

csDarmanis = sapply(as.data.frame(Z),function(z) CoreAlg(X,(z - mean(z))/sd(z))$w)
rownames(csDarmanis) = colnames(X)


#### Plotting some results ####

# Results NNLS
lsZhang_transpose <- t(lsZhang)
par(mfrow = c(2,1))
barplot(lsZhang_transpose[,1:52], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "ALS", cex.names = 0.75)
barplot(lsZhang_transpose[,52:102], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "Control", cex.names = 0.75)

# Results CIBERSORT (examples)
par(mfrow = c(2,1))
barplot(csDarmanis[,1:52], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "ALS")
barplot(csDarmanis[,52:102], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "Control")




#### To ignore ####
# # Read in IHC proportions
# 
# IHCneuro = as.matrix(read.delim('IHC.neuro.txt', sep = "\t", check.names = FALSE))[1,]
# IHCneuro = IHCneuro[intersect(names(IHCneuro),colnames(geneResid))]
# IHCastro = as.matrix(read.delim('IHC.astro.txt', sep = "\t", check.names = FALSE))[1,]
# IHCastro = IHCastro[intersect(names(IHCastro),colnames(geneResid))]
# IHCmicroglia = as.matrix(read.delim('IHC.microglia.txt', sep = "\t", check.names = FALSE))[1,]
# IHCmicroglia = IHCmicroglia[intersect(names(IHCmicroglia),colnames(geneResid))]
# IHColigo = as.matrix(read.delim('IHC.oligo.txt', sep = "\t", check.names = FALSE))[1,]
# IHColigo = IHColigo[intersect(names(IHColigo),colnames(geneResid))]
# IHCendo = as.matrix(read.delim('IHC.endo.txt', sep = "\t", check.names = FALSE))[1,]
# IHCendo = IHCendo[intersect(names(IHCendo),colnames(geneResid))]
# 
# int = unique(c(names(IHCneuro),names(IHCmicroglia),names(IHCastro),names(IHColigo),names(IHCendo)))
# IHCprops = data.frame(neuro = IHCneuro[int], astro = IHCastro[int], microglia = IHCmicroglia[int], oligo = IHColigo[int], endo = IHCendo[int], row.names = int)
