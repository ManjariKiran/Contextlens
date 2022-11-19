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
install.packages("namespace")

library(strawr)
library(devtools)
library(mariner)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)
library(InteractionSet)
library(DESeq2)

install_github('andreacirilloac/updateR')
library(updateR)
updateR()

## List loop files
loopFiles <- list.files(path = "~/contextlens/inst/extdata",pattern="Loops.txt$", full.names = TRUE)
hicFiles <- list.files(path="~contextlens/inst/extdata", pattern=".hic$", full.names = TRUE)
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
giList
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
  inner_val <- which(M >= inner) #Extract MH distance same as inner from M matrix
  inner_d <- which(!is.na(d)) #Extract MH distance same as inner from d matrix
  notNA <- intersect(inner_d,inner_val)
  new <- loop[notNA]
  return(new)
}

## Using blockApply

countMatrix_obs <-
  assay(loopCounts) |>
  aperm(c(3,4,1,2))

countMatrices <-
  assay(loopCountsoe) |>
  aperm(c(3,4,1,2))

spacings <- dim(countMatrices)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrices)[3] / nBlocks)


grid <- RegularArrayGrid(dim(countMatrices), spacings)

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

##Calculating medain of three radial distances
#X <- list(ans[[1]],ans[[2]],ans[[3]])
X <- list(ans[[1]],ans[[2]],ans[[3]])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
normalized <- apply(Y, c(1,2), median, na.rm = TRUE)

observed <- as.matrix(countMatrix_obs[11,11,,])
dim(normalized)
dim(observed)
head(loops)
## Isolate count matrix
hicFiles



norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE)
norm_MH <- norm_MH[-to_remove,]
observed <- observed[-to_remove,]
loops <- loops[-to_remove,]
colnames(observed) <- c("WT_1_1","WT_1_2","WT_2_1","WT_2_2","FS_1_1","FS_1_2","FS_2_1","FS_2_2")
colnames(normalized) <- c("WT_1_1","WT_1_2","WT_2_1","WT_2_2","FS_1_1","FS_1_2","FS_2_1","FS_2_2")
head(observed)
head(normalized)

colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep", "techrep"))
colData

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

## Run DEseq analysis
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")
summary(res1)
plotMA(res1)

## Attach DESeq2 results
mcols(loops) <- cbind(mcols(loops), observed, res1)
head(loops)
#write.table(loops,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")


# add normalization values

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

normalizationFactors(dds) <- as.matrix(norm_MH)

# perform enrichments
#dds <- DESeq(dds,betaPrior=FALSE, fitType="local")
#dds <- DESeq(dds)

## Run DEseq analysis
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_WT_vs_FS", type="apeglm")
summary(res1)
plotMA(res1)
?plotMA
plotMA(res1, ylim = NULL,
        colNonSig = "gray32", colSig = "red3", colLine = "#ff000080",
        log = "x", cex=0.45, xlab="mean expression", ylab="log fold change")
plotMA(res1, alpha = 0.1, main = "", xlab = "mean of normalized counts")

mcols(loops) <- cbind(mcols(loops), res1)
head(loops)
#write.table(loops,"WT_vs_FS_no_norm",quote=FALSE,sep="\t")
write.table(loops,"WT_vs_FS_8-10_observed",quote=FALSE,sep="\t")

l1 <- as.matrix(countMatrix_obs[,,8590,1])
l2 <- as.matrix(countMatrix_obs[,,8590,2])
l3 <- as.matrix(countMatrix_obs[,,8590,3])
l4 <- as.matrix(countMatrix_obs[,,8590,4])
l5 <- as.matrix(countMatrix_obs[,,8590,5])
l6 <- as.matrix(countMatrix_obs[,,8590,6])
l7 <- as.matrix(countMatrix_obs[,,8590,7])
l8 <- as.matrix(countMatrix_obs[,,8590,8])


df <- as.data.frame(swapAnchors(loops))[,c(1,2,8)]
tads <- GRanges(seqnames=df[[1]],
                ranges = IRanges(start=df[[2]],
                                 end = df[[3]]))

tads
WT <- Reduce("+",list(l1,l2,l3,l4))
FS <- Reduce("+",list(l5,l6,l7,l8))
dimnames(WT) <- list(seq(-100000,100000,10000),
                     seq(-100000,100000,10000))
dimnames(FS) <- list(seq(-100000,100000,10000),
                     seq(-100000,100000,10000))
WT[11,11]
ylimit <- (max(WT[11,11],FS[11,11]))*2
long_FS <-
  FS |>
  as.table() |>
  as.data.frame() |>
  setNames(c('rows','cols','counts'))
library(ggplot2)
ggplot(data = long_FS,
       mapping = aes(x= rows, y = cols, fill = counts)) +
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,ylimit) +
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
  scale_fill_distiller(palette = 'YlGnBu', direction = 1, limits =c(0,500)) +
  geom_tile() +
  theme(aspect.ratio=1, axis.text.x = element_text(angle = 45, hjust = 1))

WT
