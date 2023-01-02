# contextlens


## Overview

<img align="right" src="contextlens.png" width="150">


`contextlens` is differential loop calling package in `R`. The package identifies differential loops based on local background from Hi-C data. This package includes functions to select desired local background and utilizing those to normalize the center loop pixel to call differential loops.

## Citations
To cite `contextlens` in publications use:

`contextlens` : differential loop calling based on local context in R

## Installation

`contextlens` can be installed from `Bioconductor version 3.16`

`(R version 4.2)` as follows:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.16")

BiocManager::install("contextlens")
Example datasets and files are included with the package contextlensData:

BiocManager::install("contextlens")
```
## Usage

```R
## Load libraries and datasets
if (!requireNamespace("remotes", quietly = TRUE))
install.packages("remotes")
remotes::install_github("EricSDavis/straw/R",force =TRUE)
remotes::install_github("EricSDavis/mariner@dev", force=TRUE)
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
## Isolate count matrix

norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE) #Find the rows with NA values
low_counts <- which((rowMedians(observed) < 5) == TRUE) #Find the rows with < 5 median value
to_remove_total <- union(to_remove,low_counts)

norm_MH <- norm_MH[-to_remove_total,]
observed <- observed[-to_remove_total,]
normalized <- normalized[-to_remove_total,]
new <- countMatrix_obs[,,-to_remove_total,]
loops <- loops[-to_remove_total,]

colnames(observed) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")
colnames(normalized) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")
colnames(norm_MH) <- c("WT_0_1_1","WT_4320_1_1","WT_0_2_1","WT_4320_2_1","WT_0_1_2","WT_4320_1_2","WT_0_2_2","WT_4320_2_2","WT_0_1_3","WT_4320_1_3","WT_0_2_3","WT_4320_2_3")

colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("none","condition", "biorep", "techrep"))
colData

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
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_4320_vs_0", type="apeglm")
summary(res1)

plotMA(res1,ylim=c(-2,2),main='8-10 radial normalization',
       colSig = "skyblue",alpha = 0.05,cex.axis=1.5,cex=0.8,cex.lab=1.5)

mcols(loops) <- cbind(mcols(loops), res1)
write.table(loops,"4320_vs_0_4-6_norm",quote=FALSE,sep="\t")

```
