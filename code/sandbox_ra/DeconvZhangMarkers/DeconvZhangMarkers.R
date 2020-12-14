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

# Create a new matrix with unique genes by removing the duplicated reads (here line index) with lower expression
data_unique <- data[-c(8382, 5347, 21408, 21409, 21406, 3914, 21412, 21412, 21402, 21403, 21410, 20201, 21399,10000,21404,
                       21463, 5918, 21407, 3036, 8486, 21398, 7989, 21400, 2706, 21401, 21405, 6682, 7324, 21461, 1982,
                       21051, 21462, 21411, 3706),]
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


#### Read second marker gene ####

# Load Zhang dataset
Zhang = read.delim('Ben_Barres_HumanOnly.txt')
Zhanggenes = as.character(Zhang[,1])
Zhang = sapply(Zhang[,-1],function(x)as.numeric(as.character(x)))
rownames(Zhang) = Zhanggenes
pureSamplesZhang = list(neuro = grep('Neuro',colnames(Zhang)), astro = grep('mature.astro',colnames(Zhang)), microglia = grep('Micro',colnames(Zhang)),oligo = grep('Oligo',colnames(Zhang)),endo = grep('Endo',colnames(Zhang)))

# Reorder samples
Zhang = do.call('cbind',lapply(pureSamplesZhang,function(x)Zhang[,x]))
int = intersect(rownames(geneExpr),rownames(Zhang))
y <- DGEList(counts=Zhang[int,])
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Zhang[,1]/sum(Zhang[,1])),0.75))
Zhang <- cpm(y, log=FALSE)
pureSamplesZhang = list(neuro = grep('neuro',colnames(Zhang)), astro = grep('mature.astro',colnames(Zhang)), microglia = grep('Micro',colnames(Zhang)),oligo = grep('Oligo',colnames(Zhang)),endo = grep('Endo',colnames(Zhang)))

ZhangMean = do.call('cbind',lapply(pureSamplesZhang,function(x)if(length(x)==1){Zhang[,x]} else{rowMeans(Zhang[,x])}))
int = intersect(rownames(geneExpr),rownames(ZhangMean))
design = model.matrix(~.-1,data.frame(cell= rep(names(pureSamplesZhang),unlist(lapply(pureSamplesZhang,length)))))
colnames(design) = gsub('cell','',colnames(design))

v <- voom(Zhang[int,], design, plot=FALSE)
fit <- lmFit(v, design)
x = sapply(colnames(design),function(x)paste(x,'-(',paste(colnames(design)[colnames(design)!=x],collapse = "+"),')/4',sep = ''))
con = makeContrasts(contrasts = x,levels = colnames(design))
fit = contrasts.fit(fit,con)
fit <- eBayes(fit, robust=TRUE)
ZhangMarkers = NULL
for(i in 1:ncol(design)){
  x = fit$coef[p.adjust(fit$p.v[,i],'bonferroni')<0.05,i]
  ZhangMarkers[[colnames(design)[i]]] = names(head(sort(-x[x>0]),100))
}
unlist(lapply(ZhangMarkers,length))
ZhangMarkersFull = NULL
for(i in 1:ncol(design)){
  x = fit$coef[p.adjust(fit$p.v[,i],'fdr')<0.05,i]
  ZhangMarkersFull[[colnames(design)[i]]] =  names(sort(-x[x>0]))
}
unlist(lapply(ZhangMarkersFull,length))


#### Deconvolution ####
sce = ZhangMean
ge = geneExpr

commongenes <- intersect (rownames(ge), rownames(sce))
ge <- ge[pmatch(commongenes, rownames(ge)), ]
sce <- sce[pmatch(commongenes, rownames(sce)), ]

y <- cbind(sce, ge) # Normally it's 2^ge but it leads to 'inf' values, so maybe something is missing in the preprocess part
y <- DGEList(counts= y)
y <- calcNormFactors(y,method = 'TMM', Acutoff =  quantile(log2(Zhang[,1]/sum(Zhang[,1])),0.75))
y <- cpm(y, log=FALSE)

pure_samples = as.list(1:5)
names(pure_samples) = colnames(y)[1:5]
markers = lapply(ZhangMarkers,intersect,rownames(y))

# dtangle
dtZhang <- dtangle(log2(t(y)), pure_samples=pure_samples, markers = markers[names(pure_samples)] )$estimates

# DSA
#dsaZhang<-getcompM(geneResid,method="DSA",marker=ZhangMarkers) #unknown function
#dsaZhang <- Deconvolution(y[,6:ncol(y)], ? ) # that's the function from the package but need a weight matrix

# NNLS
dat = t(y[unlist(markers),])
lsZhang = apply(dat[-c(1:5),],1,function(x)coef(glmnet(t(dat[1:5,unlist(markers)]),x,lambda = 0, lower.limits = 0,intercept = FALSE,standardize = FALSE))[-1,])
lsZhang = t(lsZhang)/colSums(lsZhang)

# CIBERSORT
Z = y[unlist(ZhangMarkers),-c(1:5)]
X = y[unlist(ZhangMarkers),1:5]
X <- (X - mean(X))/sd(as.vector(X))
csZhang = sapply(as.data.frame(Z),function(z) CoreAlg(X,(z - mean(z))/sd(z))$w)
rownames(csZhang) = colnames(X)


#### Plotting some results (examples) ####

# Results NNLS
lsZhang_transpose <- t(lsZhang)
par(mfrow = c(2,1))
barplot(lsZhang_transpose[,1:52], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "ALS", cex.names = 0.75)
barplot(lsZhang_transpose[,52:102], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "Control", cex.names = 0.75)
# legend(0.1,0.5,c("endo","oligo","microglia","astro","neuro"),
#        cex=.8,col=c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"))

# par(mfrow = c(3,2))
# hist(lsZhang[,1])
# hist(lsZhang[,2])
# hist(lsZhang[,3])
# hist(lsZhang[,4])
# hist(lsZhang[,5])

# Results CIBERSORT
par(mfrow = c(2,1))
barplot(csZhang[,1:52], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "ALS", cex.names = 0.75)
barplot(csZhang[,52:102], col = c("lightblue", "lightcyan", "lavender", "mistyrose",  "cornsilk"),
        legend.text = TRUE, main = "Control", cex.names = 0.75)



