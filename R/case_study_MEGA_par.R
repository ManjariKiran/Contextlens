#Mega Case Study
library(strawr)
library(devtools)
library(mariner)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)
library(InteractionSet)
library(DESeq2)

## List loop files
loopFiles <- list.files(path = "~/contextlens/MEGA",pattern=".txt$", full.names = TRUE)
hicFiles <- list.files(path="~/contextlens/MEGA", pattern=".hic$", full.names = TRUE)
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
head(loops)
## Modify seqlevels of loops to match hicFiles
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Expand pixels to matrices and extract
loopCounts <-
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles,
                  norm=norm,
                  matrix=matrix)

## Do the same for observed count Expand pixels to matrices and extract
matrix <- "oe"
loopCountsoe <-
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles,
                  norm=norm,
                  matrix=matrix)

## The functions used to retreive median of signal for 1-10 radial distance
dist_diag <- function(buffer, loop, res, exclude){
  res <- 10000
  p <- pairdist(loop)
  bin <- p/res
  res1 <- res/1000
  start <- (bin)*res1
  m=(buffer*2)+1
  M <- matrix(data=NA,nrow=m,ncol=m)
  for(k in 1:m){
    rem = k-1
    for(i in k:m){
      j = i-rem
      M[i,j] <- start
    }
    start <- start - res1
  }
  start <- (bin)*res1
  for(k in 1:m){
    rem = k-1
    for(j in k:m){
      i = j-rem
      M[i,j] <- start
    }
    start <- start + res1
  }
  #  M[is.na(M)] <- 1
  M[M < exclude] <- NA
  return(M)
}
mh_index <- function(buffer, loop, inner, exclude){
  m=(buffer*2)+1
  center <- buffer+1
  M <- matrix(data=NA,nrow=m,ncol=m)
  for(j in 1:m){
    l=buffer+j
    for(i in 1:m){
      k=(m+1)-j
      if((i <= (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-k-i
      }
      if((i <= (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-j-i
      }
      if((i > (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-i-j
      }
      if((i > (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-i-k
      }
    }
  }
  d <- dist_diag(buffer, loops[i], res, exclude)
  inner_val <- which(M == inner) #Extract MH distance same as inner from M matrix
  inner_d <- which(!is.na(d)) #Extract MH distance same as inner from d matrix
  notNA <- intersect(inner_d,inner_val)
  new <- loop[notNA]
  return(new)
}

## Using blockApply

countMatrix_obs <-
  assay(loopCounts) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs)
countMatrices <-
  assay(loopCountsoe) |>
  aperm(c(3,4,1,2))

spacings <- dim(countMatrix_obs)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs)[3] / nBlocks)


grid <- RegularArrayGrid(dim(countMatrix_obs), spacings)

ans <-
  lapply(1:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs,
                 FUN = \(x) apply(x,c(3,4),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })

##Calculating median of three radial distances
X <- list(ans[[8]],ans[[9]],ans[[10]])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
normalized <- apply(Y, c(1,2), median, na.rm = TRUE)

observed <- as.matrix(countMatrix_obs[11,11,,])
dim(normalized)
dim(observed)
head(loops)
length(loops)
## Isolate count matrix

norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE)
low_counts <- which((rowMedians(observed) < 5) == TRUE)
to_remove_total <- union(to_remove,low_counts)

norm_MH <- norm_MH[-to_remove_total,]
observed <- observed[-to_remove_total,]
normalized <- normalized[-to_remove_total,]
new <- countMatrix_obs[,,-to_remove_total,]
dim(norm_MH)
dim(observed)
#Remove Loops with a median count of 5 counts or less

dim(observed)
loops <- loops[-to_remove_total,]
colnames(observed) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")
colnames(normalized) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")
colnames(norm_MH) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")
head(observed)
head(normalized)
head(norm_MH)

colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("none","condition", "biorep", "techrep"))
colData

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition + biorep)

default <- counts(dds,normalized = FALSE)
head(default)
default_sizefactor <- counts(dds, normalized =TRUE)
sizeFactors(dds)
head(default_sizefactor)
## Run DEseq analysis

dds <- DESeq(dds)
resultsNames(dds)
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_4320_vs_0", type="apeglm")
summary(res1)
res2 <- na.omit(res1)
res3 <- res2[((res2$padj < 0.05) == "TRUE"),]
res4 <- res3[((abs(res3$log2FoldChange) > 0.585) == "TRUE"),]
dim(res4)
plotMA(res1,ylim=c(-2,2),main='No normalization',
       colSig = "skyblue",alpha = 0.05,cex.axis=1.5,cex=0.8,cex.lab=1.5)
?geneplotter::plotMA
dim(res1)
## Attach DESeq2 results
mcols(loops) <- cbind(mcols(loops), res1)
head(loops)
write.table(loops,"4320_vs_0_no_norm",quote=FALSE,sep="\t")


# add normalization values

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition + biorep)

normalizationFactors(dds) <- as.matrix(norm_MH)
default <- counts(dds,normalized = FALSE)
head(default)
default_normfactor <- counts(dds, normalized =TRUE)
head(default_normfactor)

# perform enrichments
#dds <- DESeq(dds,betaPrior=FALSE, fitType="local")
dds <- DESeq(dds)
normalizationFactors(dds) <- as.matrix(norm_MH)

## Run DEseq analysis
res2 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_4320_vs_0", type="apeglm")
summary(res1)
res2 <- na.omit(res1)
res3 <- res2[((res2$padj < 0.05) == "TRUE"),]
res4 <- res3[((abs(res3$log2FoldChange) > 0.585) == "TRUE"),]
dim(res4)
plotMA(res1,ylim=c(-2,2),main='8-10 radial normalization',
       colSig = "skyblue",alpha = 0.05,cex.axis=1.5,cex=0.8,cex.lab=1.5)

head(res1)
length(loops)
mcols(loops) <- cbind(mcols(loops), res1,res2)
head(loops)
#write.table(loops,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")
write.table(loops,"4320_vs_0_8-10_norm",quote=FALSE,sep="\t")

loop <- 1515
observed[loop,]
normalized[loop,]
norm_MH[loop,]
default_normfactor[loop,]
dim(observed)
dim(normalized)
dim(norm_MH)
dim(default_normfactor)
dim(new)
hicFiles
l1 <- as.matrix(new[,,loop,1])/norm_MH[loop,1]
l2 <- as.matrix(new[,,loop,3])/norm_MH[loop,3]
l3 <- as.matrix(new[,,loop,5])/norm_MH[loop,5]
l4 <- as.matrix(new[,,loop,7])/norm_MH[loop,7]
l5 <- as.matrix(new[,,loop,9])/norm_MH[loop,9]
l6 <- as.matrix(new[,,loop,11])/norm_MH[loop,11]
l7 <- as.matrix(new[,,loop,2])/norm_MH[loop,2]
l8 <- as.matrix(new[,,loop,4])/norm_MH[loop,4]
l9 <- as.matrix(new[,,loop,6])/norm_MH[loop,6]
l10 <- as.matrix(new[,,loop,8])/norm_MH[loop,8]
l11 <- as.matrix(new[,,loop,10])/norm_MH[loop,10]
l12 <- as.matrix(new[,,loop,12])/norm_MH[loop,12]

l1 <- as.matrix(new[,,loop,1])
l2 <- as.matrix(new[,,loop,3])
l3 <- as.matrix(new[,,loop,5])
l4 <- as.matrix(new[,,loop,7])
l5 <- as.matrix(new[,,loop,9])
l6 <- as.matrix(new[,,loop,11])
l7 <- as.matrix(new[,,loop,2])
l8 <- as.matrix(new[,,loop,4])
l9 <- as.matrix(new[,,loop,6])
l10 <- as.matrix(new[,,loop,8])
l11 <- as.matrix(new[,,loop,10])
l12 <- as.matrix(new[,,loop,12])


WT_0 <- Reduce("+",list(l1,l2,l3,l4,l5,l6))
WT_72 <- Reduce("+",list(l7,l8,l9,l10,l11,l12))
dimnames(WT_0) <- list(seq(-100000,100000,10000),
                       seq(-100000,100000,10000))
dimnames(WT_72) <- list(seq(-100000,100000,10000),
                        seq(-100000,100000,10000))
#ylimit <- (max(WT[11,11],FS[11,11]))
ylimit <- (max(c(WT_0,WT_72),na.rm=TRUE))
#ylimit <- ylimit*2
ylimit
long_WT_0 <-
  WT_0 |>
  as.table() |>
  as.data.frame() |>
  setNames(c('rows','cols','counts'))
library(ggplot2)
ggplot(data = long_WT_0,
       mapping = aes(x= rows, y = cols, fill = counts)) +
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,ylimit)) +
  geom_tile() +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1))
long_WT_72 <-
  WT_72 |>
  as.table() |>
  as.data.frame() |>
  setNames(c('rows','cols','counts'))
library(ggplot2)
ggplot(data = long_WT_72,
       mapping = aes(x= rows, y = cols, fill = counts)) +
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,ylimit)) +
  geom_tile() +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1))
