##Program to simulate loops from given Hi-C files
library(strawr)
library(devtools)
library(mariner)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)
library(InteractionSet)
library(DESeq2)

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

##Set a distance threshold
loops10kb <- loops[which(pairdist(loops) > 100000),]

## Expand pixels to matrices and extract
loopCounts <-
  pixelsToMatrices(x=loops10kb, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5:8],
                  norm=norm,
                  matrix=matrix)

## Using blockApply

countMatrix_obs <-
  assay(loopCounts) |>
  aperm(c(3,4,1,2))

## Expand pixels to matrices and generate APA
loopCount_APA <-
  pixelsToMatrices(x=loops10kb, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5:8],
                  norm=norm,
                  matrix=matrix,
                  onDisk = F) |>
                  aggHicMatrices(FUN=sum) |>
                  as.matrix()

dim(loopCount_APA)

#Calculate decay for each pixel from  center

m=(buffer*2)+1
center <- buffer+1
decay_matrix <- matrix(data=NA,nrow=m,ncol=m)
for(i in 1:21){
  for(j in 1:21){
    decay_matrix[i,j] <- loopCount_APA[i,j]/loopCount_APA[center,center]
  }
}

decay_matrix_subs <- (1-decay_matrix)
decay_matrix_subs


##Calculate mean and sd of center pixel
observed <- as.matrix(countMatrix_obs[11,11,,])
colnames(observed) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2")
dim(observed)

loop_mean <- rowMeans(observed)
loop_sd <- rowVars(observed)
observed_stat <- cbind(observed,round(loop_mean))
colnames(observed_stat) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","mean")
number <- sort(unique(observed_stat[,5]))
ref <- data.frame("mean","sd")
colnames(ref) <- c("mean","sd")
head(ref)
for(i in 2:length(number)){
  x <- mean(observed_stat[which(observed_stat[,5] == number[i]),],na.rm=T)
  y <- sd(observed_stat[which(observed_stat[,5] == number[i]),],na.rm=T)
  ref[i,] <- cbind(number[i],round(y,3))
}

##Print random numbers for given mean and std
center_df <- data.frame("rep1","rep2","rep3","rep4","rep5","rep6")
colnames(center_df) <- data.frame("rep1","rep2","rep3","rep4","rep5","rep6")
for(i in 1:10000){
   random <-  as.numeric(sample(number, 1, replace = TRUE))
   n_sd <- as.numeric(ref[which(ref$mean == random),]$sd)
   center_df[i,] <- c(round(rnorm(n=6,mean = random, sd= n_sd),0))
}
dat <- as.data.frame(sapply(center_df, as.numeric))
str(dat)
WT_mean <- round(rowMeans(dat[,1:3]))
center_df <- cbind(center_df,WT_mean)
head(center_df)
##Select random loops for fold change
head(center_df)
center_df2 <- center_df[,4:6]
head(center_df2)
##Create 10% fold change loops
de <- sample(1:10000,size = 400, replace=F)
#de_half <- sample(de,size=400,replace=F)
#de_rem <- de[!de %in% de_half]
#de_two <- sample(de_rem,size=300,replace=F)
#de_rem <- de_rem[!de_rem %in% de_two]
#de_three <- sample(de_rem,size=200,replace=F)
de_four_cnv <- sample(de,size = 200, replace =F)
de_rem <- de[!de %in% de_four_cnv]
#de_rem <- de_rem[!de_rem %in% de_three]
de_four <- de_rem

write.csv(de_four_cnv,"de_four_cnv")
write.csv(de_four,"de_four")

##Set these as center pixel and create local background
Rep1 <- vector("list",length = 10000)
Rep2 <- vector("list",length = 10000)
Rep3 <- vector("list",length = 10000)
Rep4 <- vector("list",length = 10000)
Rep5 <- vector("list",length = 10000)
Rep6 <- vector("list",length = 10000)
Rep7 <- vector("list",length = 10000)
m=(buffer*2)+1
for(i in 1:10000){
  M1 <- matrix(data=NA,nrow=m,ncol=m)
  M2 <- matrix(data=NA,nrow=m,ncol=m)
  M3 <- matrix(data=NA,nrow=m,ncol=m)
  M4 <- matrix(data=NA,nrow=m,ncol=m)
  M5 <- matrix(data=NA,nrow=m,ncol=m)
  M6 <- matrix(data=NA,nrow=m,ncol=m)
  M7 <- matrix(data=NA,nrow=m,ncol=m)
      M1 <- round(as.numeric(center_df[i,1]) - decay_matrix_subs*as.numeric(center_df[i,1]))
      M2 <- round(as.numeric(center_df[i,2]) - decay_matrix_subs*as.numeric(center_df[i,2]))
      M3 <- round(as.numeric(center_df[i,3]) - decay_matrix_subs*as.numeric(center_df[i,3]))
      M4 <- round(as.numeric(center_df[i,4]) - decay_matrix_subs*as.numeric(center_df[i,4]))
      M5 <- round(as.numeric(center_df[i,5]) - decay_matrix_subs*as.numeric(center_df[i,5]))
      M6 <- round(as.numeric(center_df[i,6]) - decay_matrix_subs*as.numeric(center_df[i,6]))
      M7 <- round(as.numeric(center_df[i,7]) - decay_matrix_subs*as.numeric(center_df[i,7]))
Rep1[[i]] <- M1
Rep2[[i]] <- M2
Rep3[[i]] <- M3
Rep4[[i]] <- M4
Rep5[[i]] <- M5
Rep6[[i]] <- M6
Rep7[[i]] <- M7
}
ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)
length(de_four)
de_four
Rep1[[6923]]
Rep4[[6923]]
for(i in 1:length(de_four)){
  x <- de_four[i]
  sides <- median(c(Rep7[[x]][1,],Rep7[[x]][,1],Rep7[[x]][m,],Rep7[[x]][,m]),na.rm=TRUE)
  new_center4 <- 4*as.numeric(center_df[x,4])
  new_center5 <- 4*as.numeric(center_df[x,5])
  new_center6 <- 4*as.numeric(center_df[x,6])
  alpha4 = (new_center4 - sides)/(as.numeric(center_df[x,7]) - sides)
  alpha5 = (new_center5 - sides)/(as.numeric(center_df[x,7]) - sides)
  alpha6 = (new_center6 - sides)/(as.numeric(center_df[x,7]) - sides)
  Rep4[[x]] <- round(new_center4 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha4))
  Rep5[[x]] <- round(new_center5 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha5))
  Rep6[[x]] <- round(new_center6 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha6))
}

ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)

for(i in 1:length(de_four_cnv)){
  x <- de_four_cnv[i]
  new_center4 <- 4*as.numeric(center_df[x,4])
  new_center5 <- 4*as.numeric(center_df[x,5])
  new_center6 <- 4*as.numeric(center_df[x,6])
  Rep4[[x]] <- round(new_center4 - decay_matrix_subs*new_center4)
  Rep5[[x]] <- round(new_center5 - decay_matrix_subs*new_center5)
  Rep6[[x]] <- round(new_center6 - decay_matrix_subs*new_center6)
}

ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)

ALL_simulated <- array(unlist(ALL),dim=c(21,21,10000,6))
ALL <- DelayedArray::DelayedArray(ALL_simulated)

spacings <- dim(ALL)
nBlocks <- 20
spacings[3] <- ceiling(dim(ALL)[3] / nBlocks)

grid <- RegularArrayGrid(dim(ALL), spacings)
grid
ans <-
  lapply(8:10, \(i) {
    tmp <-
      blockApply(x = ALL,
                 FUN = \(x) apply(x,c(3,4),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })

##Calculating median of three radial distances
X <- list(ans[[1]],ans[[2]],ans[[3]])
Y <- do.call(cbind, X)
Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
normalized <- apply(Y, c(1,2), median, na.rm = TRUE)
dim(ALL)

observed <- ALL[11,11,,]
head(observed)

## Isolate count matrix
loop <- data.frame(seq(1:10000))
loop1 <- cbind(loop,observed)
colnames(loop1) <- c("loop","WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE)
norm_MH <- norm_MH[-to_remove,]
observed <- observed[-to_remove,]
normalized <- normalized[-to_remove,]
new <- ALL[,,-to_remove,]
loops <- loop1[-to_remove,]
l <- subset(loops, loop %in% de)

dim(new)
dim(loops)

##Add column names
colnames(observed) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
colnames(normalized) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
colnames(norm_MH) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
head(observed)
head(normalized)
head(norm_MH)
##No normalization
colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep"))
colData

dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)
dds$condition <- factor(dds$condition, levels=c("WT","TRT"))
default <- counts(dds,normalized = FALSE)
#default_sizefactor <- counts(dds, normalized =TRUE)
sizeFactors(dds)
head(default_sizefactor)
## Run DEseq analysis

dds <- DESeq(dds)
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_TRT_vs_WT", type="apeglm")
summary(res1)
file <- cbind(loops,res1)
head(file)
write.csv(file,"no_normalization.csv")
##With normalization
dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

normalizationFactors(dds) <- as.matrix(norm_MH)
default <- counts(dds,normalized = FALSE)
head(default)
default_normfactor <- counts(dds, normalized =TRUE)
head(default_normfactor)
head(norm_MH)

# perform enrichments

dds$condition <- factor(dds$condition, levels=c("WT","TRT"))
normalizationFactors(dds) <- as.matrix(norm_MH)

## Run DEseq analysis
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_TRT_vs_WT", type="apeglm")
summary(res1)
file1 <- cbind(loops,res1)
head(file1)
write.csv(file1,"normalization.csv")
subset(file1, loop %in% de_four) |>
    subset(log2FoldChange < 2) |>
    subset(padj > 0.1)
subset(file, loop %in% de_four) |>
  subset(log2FoldChange < 2) |>
  subset(padj > 0.1)
center_df[4155,]


