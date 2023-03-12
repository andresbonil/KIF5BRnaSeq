library(gplots)
library(Rsubread)
library(edgeR)
library(RColorBrewer)
library(Glimma)
library(biomaRt)
wd <- getwd()

#7709
data <- read.delim("7709_counts.txt", row.names = NULL, stringsAsFactors = FALSE)

#Find external gene names from gene_id
geneids <- data.frame(data$GeneID)
ensembl <- useEnsembl(biomart = "genes")
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attr <- listAttributes(ensembl.con)

genesymbols <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                     filters = "ensembl_gene_id",
                     values = geneids$data.GeneID,
                     mart = ensembl.con)

sortdata <- data[order(data$GeneID),]
sortnames <- genesymbols[order(genesymbols$ensembl_gene_id),]
sortdata <- cbind(sortdata, sortnames$external_gene_name)

sortdata$`sortnames$external_gene_name` <- make.unique(sortdata$`sortnames$external_gene_name`)
rownames(sortdata) <- sortdata$`sortnames$external_gene_name`
count.data <- sortdata[3:9]

#7804
data <- read.delim("7804_counts.txt", row.names = NULL, stringsAsFactors = FALSE)

#Find external gene names from gene_id
geneids <- data.frame(data$GeneID)
ensembl <- useEnsembl(biomart = "genes")
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

attr <- listAttributes(ensembl.con)

genesymbols <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                     filters = "ensembl_gene_id",
                     values = geneids$data.GeneID,
                     mart = ensembl.con)

sortdata <- data[order(data$GeneID),]
sortnames <- genesymbols[order(genesymbols$ensembl_gene_id),]
sortdata <- cbind(sortdata, sortnames$external_gene_name)

sortdata$`sortnames$external_gene_name` <- make.unique(sortdata$`sortnames$external_gene_name`)
rownames(sortdata) <- sortdata$`sortnames$external_gene_name`
count.data.2 <- sortdata[3:9]

countpool <- cbind(count.data, count.data.2)

grouppool <- c('WT', 'WT_si5B1', 'WT_si5B2', 'WT_si1B1', 'WT_si1B4', 'WT_siScramble', 'KIF5BKO', 'WT', 'WT_si5B1', 'WT_si5B2', 'WT_si1B1', 'WT_si1B4', 'WT_siScramble', 'KIF5BKO')
experimentpool <- c('one','one','one','one','one','one','one','two','two','two','two','two','two','two')
group <- c('1_1', '1_2', '1_3', '1_4', '1_5', '1_6', '1_7')
group.2 <- c('2_1', '2_2', '2_3', '2_4', '2_5', '2_6', '2_7')

dge.set <- DGEList(counts = count.data,
                   samples = group,
                   genes = sortdata[-c(2:9)])
dge.set.2 <- DGEList(counts = count.data.2,
                   samples = group.2,
                   genes = sortdata[-c(2:9)])
dge.set.3 <- DGEList(counts = countpool,
                     group = experimentpool,
                     samples = grouppool,
                     genes = sortdata[-c(2:9)])

samplenames <- paste(gsub('_','',substr(colnames(dge.set),1,8)),group,sep='_')
samplenames2 <- paste(gsub('_','',substr(colnames(dge.set.2),1,8)),group,sep='_')
samplenamespool <- paste(gsub('_','',substr(colnames(dge.set.3),1,8)),group,sep='_')

lcpm.prefilter <- cpm(dge.set,log=T)
lcpm.prefilter.2 <- cpm(dge.set.2,log=T)
lcpm.prefilter.3 <- cpm(dge.set.3, log=T)

L <- mean(dge.set$samples$lib.size) * 1e-6
M <- median(dge.set$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm.prefilter)

L <- mean(dge.set.2$samples$lib.size) * 1e-6
M <- median(dge.set.2$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm.prefilter.2)

L <- mean(dge.set.3$samples$lib.size) * 1e-6
M <- median(dge.set.3$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm.prefilter.3)


table(rowSums(dge.set$counts==0)==7)
table(rowSums(dge.set.2$counts==0)==7)
table(rowSums(dge.set.3$counts==0)==14)

keep.exprs <- filterByExpr(dge.set, group=group)
x <- dge.set[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge.set)
dim(x)

keep.exprs <- filterByExpr(dge.set.2, group=group)
y <- dge.set.2[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge.set.2)
dim(y)


keep.exprs <- filterByExpr(dge.set.3, group=grouppool)
z <- dge.set.3[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge.set.3)
dim(z)

#7709 QC figures  

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm.prefilter[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm.prefilter[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors 

#7709 MDS plot
group.index <- rep(0,length(group))
for(i in 1:length(unique(group))){
  group.index[which(group==unique(group)[i])] <- i
}
glMDSPlot(lcpm, top=NULL, groups=group, labels=samplenames, html='All_Genes', folder='MDS_Plots', launch=TRUE)
glMDSPlot(lcpm, top=500, groups=group, gene.selection='common', folder='MDS_Plots',
          labels=samplenames, html='Top500_variable_genes', launch=TRUE)

#7804 QC figures 

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm.prefilter.2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm.prefilter.2[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
x <- calcNormFactors(y, method = "TMM")
x$samples$norm.factors 

#7804 MDS plot 
group.index <- rep(0,length(group.2))
for(i in 1:length(unique(group.2))){
  group.index[which(group.2==unique(group.2)[i])] <- i
}
glMDSPlot(lcpm, top=NULL, groups=group.2, labels=samplenames, html='All_Genes', folder='MDS_Plots', launch=TRUE)
glMDSPlot(lcpm, top=500, groups=group.2, gene.selection='common', folder='MDS_Plots',
          labels=samplenames, html='Top500_variable_genes', launch=TRUE)


#Pooled QC figures 

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(z)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm.prefilter.3[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm.prefilter.3[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

lcpm <- cpm(z, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
z <- calcNormFactors(z, method = "TMM")
z$samples$norm.factors 

#Pooled MDS plot 
group.index <- rep(0,length(grouppool))
for(i in 1:length(unique(grouppool))){
  group.index[which(group.2==unique(grouppool)[i])] <- i
}
glMDSPlot(lcpm, top=NULL, groups=grouppool, labels=samplenamespool, html='All_Genes', folder='MDS_Plots', launch=TRUE)
glMDSPlot(lcpm, top=500, groups=grouppool, gene.selection='common', folder='MDS_Plots',
          labels=samplenamespool, html='Top500_variable_genes', launch=TRUE)

#Differential Expression 

design <- model.matrix(~0+group)
design2 <- model.matrix(~0+group.2)
designpool <- model.matrix(~0+grouppool+experimentpool)
colnames(design) <- gsub("group1_", "f", colnames(design))
colnames(design2) <- gsub("group.22_" ,"s", colnames(design2))
colnames (designpool) <- gsub("grouppool", "", colnames(designpool))
contr.matrix <- makeContrasts(
  WTvsWT_si5B1 = WT - WT_si5B1,
  WTvsWT_si5B2 = WT - WT_si5B2,
  WTvsWT_si1B1 = WT - WT_si1B1,
  WTvsWT_si1B4 = WT - WT_si1B4,
  WTvsWT_siScramble = WT - WT_siScramble,
  WTvsKIF5BKO = WT - KIF5BKO,
  WT_si5B1vsWT_si5B2 = WT_si5B1 - WT_si5B2,
  WT_si5B1vsWT_si1B1 = WT_si5B1 - WT_si1B1,
  WT_si5B1vsWT_si1B4 = WT_si5B1 - WT_si1B4,
  WT_si5B1vsWT_siScramble = WT_si5B1 - WT_siScramble,
  WT_si5B1vsKIF5BKO = WT_si5B1 - KIF5BKO,
  WT_si5B2vsWT_si1B1 = WT_si5B2 - WT_si1B1,
  WT_si5B2vsWT_si1B4 = WT_si5B2 - WT_si1B4,
  WT_si5B2vsWT_siScramble = WT_si5B2 - WT_siScramble,
  WT_si5B2vsKIF5BKO = WT_si5B2 - KIF5BKO,
  WT_si1B1vsWT_si1B4 = WT_si1B1 - WT_si1B4,
  WT_si1B1vsWT_siScramble = WT_si1B1 - WT_siScramble,
  WT_si1B1vsKIF5BKO = WT_si1B1 - KIF5BKO,
  WT_si1B4vsWT_siScramble = WT_si1B4 - WT_siScramble,
  WT_si1B4vsKIF5BKO = WT_si1B4 - KIF5BKO,
  WT_siScramblevsKIF5BKO = WT_siScramble - KIF5BKO,
  levels = colnames(designpool))

#Voom plot, z is pooled data
v <- voom(z, designpool, plot=TRUE)

vfit <- lmFit(v, designpool)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=2, status=dt, main=colnames(tfit)[2], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=3, status=dt, main=colnames(tfit)[3], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=4, status=dt, main=colnames(tfit)[4], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=5, status=dt, main=colnames(tfit)[5], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=6, status=dt, main=colnames(tfit)[6], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=7, status=dt, main=colnames(tfit)[7], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=8, status=dt, main=colnames(tfit)[8], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=9, status=dt, main=colnames(tfit)[9], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=10, status=dt, main=colnames(tfit)[10], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=11, status=dt, main=colnames(tfit)[11], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=12, status=dt, main=colnames(tfit)[12], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=13, status=dt, main=colnames(tfit)[13], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=14, status=dt, main=colnames(tfit)[14], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=15, status=dt, main=colnames(tfit)[15], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=16, status=dt, main=colnames(tfit)[16], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=17, status=dt, main=colnames(tfit)[17], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=18, status=dt, main=colnames(tfit)[18], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=19, status=dt, main=colnames(tfit)[19], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=20, status=dt, main=colnames(tfit)[20], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)
glMDPlot(tfit, coef=21, status=dt, main=colnames(tfit)[21], side.main="GeneName", counts=lcpm, groups=grouppool, launch=TRUE)

