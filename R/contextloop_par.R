if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("EricSDavis/straw/R",force =TRUE)
remotes::install_github("EricSDavis/mariner@dev", force=TRUE)
remotes::install_github("EricSDavis/mariner@dev")
install.packages("SummarizedExperiment")
install.packages('devtools')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
#install.packages("namespace")

library(strawr)
library(devtools)
library(mariner)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)
library(InteractionSet)
library(DESeq2)

#install_github('andreacirilloac/updateR')
#library(updateR)
#updateR()

## List loop files
loopFiles <- list.files(path = "~/contextlens/inst/extdata",pattern="Loops.txt$", full.names = TRUE)
hicFiles <- list.files(path="~/contextlens/inst/extdata", pattern=".hic$", full.names = TRUE)
hicFiles
loopFiles
## Define parameters
buffer <- 10
res <- 10e03
norm <- "NONE"
matrix <- "observed"

## Read in loop files as GInteractions
giList <-
  lapply(loopFiles, read.table, header=TRUE) |>
  lapply(as_ginteractions) |>
  setNames(gsub(".*//(.*)_.*", "\\1", loopFiles))
## Merge into a single set of loops
loops <- mergePairs(x = giList, radius = 20e03, column = "APScoreAvg", selectMax = TRUE)
head(loops)

## Change 5kb loops to 10kb
loops <- binPairs(x=loops, binSize=res)

## Modify seqlevels of loops to match hicFiles
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Expand pixels to matrices and extract
loopCounts <-
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles,
                  norm=norm,
                  matrix=matrix)

countMatrix_obs <-
  assay(loopCounts) |>
  aperm(c(3,4,1,2))

##Function to get median of pixel parallel to diagnol at every MH distance
radial_par <- function(buffer, loop, inner){
  center <- buffer+1
  output <- c()
  output1 <- loop[center+inner,center+inner]
  output2 <- loop[center-inner,center-inner]
  output <- c(output1,output2)
  return(output)
}

## Using blockApply
spacings <- dim(countMatrix_obs)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs)[3] / nBlocks)

grid <- RegularArrayGrid(dim(countMatrix_obs), spacings)

ans <-
  lapply(1:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs,
                 FUN = \(x) apply(x,c(3,4),FUN = \(z) {
                   median(radial_par(buffer = buffer, loop=z, inner=i))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })


##Calculating medain of three radial distances
X <- list(ans[[8]],ans[[9]],ans[[10]])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
normalized <- apply(Y, c(1,2), median, na.rm = TRUE)

loops <- mergePairs(x = giList, radius = 20e03, column = "APScoreAvg", selectMax = TRUE)
loops <- binPairs(x=loops, binSize=res)
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
observed <- as.matrix(countMatrix_obs[11,11,,])
dim(normalized)
dim(observed)
head(loops)
## Isolate count matrix

norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE)
norm_MH <- norm_MH[-to_remove,]
observed <- observed[-to_remove,]
normalized <- normalized[-to_remove,]
new <- countMatrix_obs[,,-to_remove,]
loops <- loops[-to_remove,]
colnames(observed) <- c("WT_1_1","WT_1_2","WT_2_1","WT_2_2","FS_1_1","FS_1_2","FS_2_1","FS_2_2")
colnames(normalized) <- c("WT_1_1","WT_1_2","WT_2_1","WT_2_2","FS_1_1","FS_1_2","FS_2_1","FS_2_2")
colnames(norm_MH) <- c("WT_1_1","WT_1_2","WT_2_1","WT_2_2","FS_1_1","FS_1_2","FS_2_1","FS_2_2")
#colnames(observed) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","WT_1_1","WT_1_2","WT_2_1","WT_2_2")
#colnames(normalized) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","WT_1_1","WT_1_2","WT_2_1","WT_2_2")
#colnames(norm_MH) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","WT_1_1","WT_1_2","WT_2_1","WT_2_2")
dim(observed)
dim(normalized)
dim(norm_MH)
length(loops)

colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep", "techrep"))
colData

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

#default <- counts(dds,normalized = FALSE)
#head(default)
#default_sizefactor <- counts(dds, normalized =TRUE)
#sizeFactors(dds)
#head(default_sizefactor)
## Run DEseq analysis
dds$condition
dds$condition <- factor(dds$condition, level = c("WT","FS"))
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_FS_vs_WT", type="apeglm")
summary(res1)
plotMA(res1,ylim=c(-2,2),main='No normalization',
      colSig = "skyblue",alpha = 0.1,cex.axis=1.5,cex=0.8,cex.lab=1.5,
       cex.main=1.5,ylab="Log2 Fold Change",xlab="Mean of Normalized Counts",cex.main=1.5)

## Separate WT/FS-specific loops
#res2 <- na.omit(res1)
#res3 <- res2[((res2$padj < 0.1) == "TRUE"),]
#res4 <- res3[((abs(res3$log2FoldChange) > 0.585) == "TRUE"),]
#dim(res3)
## Attach DESeq2 results
mcols(loops) <- cbind(mcols(loops), observed, res1)
head(loops)

#write.table(loops,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")


# add normalization values

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

#normalizationFactors(dds) <- as.matrix(norm_MH)
#default <- counts(dds,normalized = FALSE)
#head(default)
default_normfactor <- counts(dds, normalized =TRUE)
#head(default_normfactor)

# perform enrichments
#dds <- DESeq(dds,betaPrior=FALSE, fitType="local")
#dds <- DESeq(dds)
normalizationFactors(dds) <- as.matrix(norm_MH)
dds$condition <- factor(dds$condition, level = c("WT","FS"))
## Run DEseq analysis
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_FS_vs_WT", type="apeglm")
summary(res1)
#plotMA(res1)
#?plotMA
plotMA(res1,ylim=c(-2,2),main='8-10 normaliztion',
       colSig = "skyblue",alpha = 0.1,cex.axis=1.5,cex=0.8,cex.lab=1.5,
       cex.main=1.5,ylab="Log2 Fold Change",xlab="Mean of Normalized Counts",cex.main=1.5)

res2 <- na.omit(res1)
res3 <- res2[((res2$padj < 0.1) == "TRUE"),]
up <- res3[(((res3$log2FoldChange) >= 0) == "TRUE"),]
down <- res3[(((res3$log2FoldChange) < 0) == "TRUE"),]
dim(up)
dim(down)

mcols(loops) <- cbind(mcols(loops), res1)
head(loops)
#write.table(loops,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")
write.table(loops,"WT_vs_FS_8-10",quote=FALSE,sep="\t")

loop <- 876
observed[loop,]
normalized[loop,]
norm_MH[loop,]
default_normfactor[loop,]
new <- countMatrix_obs[,,-to_remove,]
l1 <- as.matrix(new[,,loop,1])/norm_MH[loop,1]
l2 <- as.matrix(new[,,loop,2])/norm_MH[loop,2]
l3 <- as.matrix(new[,,loop,3])/norm_MH[loop,3]
l4 <- as.matrix(new[,,loop,4])/norm_MH[loop,4]
l5 <- as.matrix(new[,,loop,5])/norm_MH[loop,5]
l6 <- as.matrix(new[,,loop,6])/norm_MH[loop,6]
l7 <- as.matrix(new[,,loop,7])/norm_MH[loop,7]
l8 <- as.matrix(new[,,loop,8])/norm_MH[loop,8]
l1 <- as.matrix(new[,,loop,1])
l2 <- as.matrix(new[,,loop,2])
l3 <- as.matrix(new[,,loop,3])
l4 <- as.matrix(new[,,loop,4])
l5 <- as.matrix(new[,,loop,5])
l6 <- as.matrix(new[,,loop,6])
l7 <- as.matrix(new[,,loop,7])
l8 <- as.matrix(new[,,loop,8])

WT <- Reduce("+",list(l1,l2,l3,l4))
FS <- Reduce("+",list(l5,l6,l7,l8))
dimnames(WT) <- list(seq(-100000,100000,10000),
                     seq(-100000,100000,10000))
dimnames(FS) <- list(seq(-100000,100000,10000),
                     seq(-100000,100000,10000))
#ylimit <- (max(WT[11,11],FS[11,11]))
ylimit <- (max(c(WT,FS),na.rm=TRUE))
#ylimit <- ylimit*2
ylimit
long_FS <-
  FS |>
  as.table() |>
  as.data.frame() |>
  setNames(c('rows','cols','counts'))
library(ggplot2)
ggplot(data = long_FS,
       mapping = aes(x= rows, y = cols, fill = counts)) +
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,ylimit)) +
  geom_tile() +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1))
long_WT <-
  WT |>
  as.table() |>
  as.data.frame() |>
  setNames(c('rows','cols','counts'))
library(ggplot2)
ggplot(data = long_WT,
       mapping = aes(x= rows, y = cols, fill = counts)) +
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,ylimit)) +
  geom_tile() +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1))
