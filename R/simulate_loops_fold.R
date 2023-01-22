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
for(i in 1:m){
  for(j in 1:m){
    decay_matrix[i,j] <- loopCount_APA[i,j]/loopCount_APA[center,center]
  }
}

decay_matrix_subs <- (1-decay_matrix)
decay_matrix_subs


##Calculate mean and sd of center pixel
observed <- as.matrix(countMatrix_obs[11,11,,])
#observed <- as.matrix(countMatrix_obs[16,16,,])
colnames(observed) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2")
observed

loop_mean <- rowMeans(observed)
loop_sd <- rowVars(observed)
observed_stat <- cbind(observed,round(loop_mean))
colnames(observed_stat) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","mean")
number <- sort(unique(observed_stat[,5]))
ref <- data.frame("mean","sd")
colnames(ref) <- c("mean","sd")
head(ref)
number <- number[-1]
number
for(i in 1:length(number)){
  x <- mean(observed_stat[which(observed_stat[,5] == number[i]),],na.rm=T)
  y <- sd(observed_stat[which(observed_stat[,5] == number[i]),],na.rm=T)
  ref[i,] <- cbind(number[i],round(y,3))
}
head(ref)
number
##Print random numbers for given mean and std
center_df <- data.frame("rep1","rep2","rep3","rep4","rep5","rep6")
colnames(center_df) <- data.frame("rep1","rep2","rep3","rep4","rep5","rep6")
for(i in 1:10000){
   random <-  as.numeric(sample(number, 1, replace = TRUE))
   n_sd <- as.numeric(ref[which(ref$mean == random),]$sd)
   center_df[i,] <- c(round(rnorm(n=6,mean = random, sd= n_sd),0))
}
dat <- as.data.frame(sapply(center_df, as.numeric))
head(dat)
WT_mean <- round(rowMeans(dat[,1:3]))
center_df <- cbind(center_df,WT_mean)
head(center_df)
##Select random loops for fold change
head(center_df)
center_df2 <- center_df[,4:6]
head(center_df2)
##Create 10% fold change loops
de <- sample(1:10000,size = 400, replace=F)
de_half <- sample(de,size=400,replace=F)
de_rem <- de[!de %in% de_half]
de_two <- sample(de_rem,size=300,replace=F)
de_rem <- de_rem[!de_rem %in% de_two]
de_three <- sample(de_rem,size=200,replace=F)
#de_rem <- de[!de %in% de_four_cnv]
de_rem <- de_rem[!de_rem %in% de_three]
de_four_cnv <- sample(de,size = 200, replace =F)
de_four <- de_rem
length(de)
length(de_half)
length(de_two)
length(de_three)
length(de_four)
de_half_cnv <- sample(de_half,size=200,replace=F)
de_half <- de_half[!de_half %in% de_half_cnv]
de_two_cnv <- sample(de_two,size=150,replace=F)
de_two <- de_two[!de_two %in% de_two_cnv]
de_three_cnv <- sample(de_three,size=100,replace=F)
de_three <- de_three[!de_three %in% de_three_cnv]
de_four_cnv <- sample(de_four,size=50,replace=F)
de_four <- de_four[!de_four %in% de_four_cnv]



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
#Rep1[[6923]]
#Rep4[[6923]]
for(i in 1:length(de_half)){
  x <- de_half[i]
  sides <- median(c(Rep7[[x]][1,],Rep7[[x]][,1],Rep7[[x]][m,],Rep7[[x]][,m]),na.rm=TRUE)
  new_center4 <- 1.5*as.numeric(center_df[x,4])
  new_center5 <- 1.5*as.numeric(center_df[x,5])
  new_center6 <- 1.5*as.numeric(center_df[x,6])
  alpha4 = (new_center4 - sides)/(as.numeric(center_df[x,7]) - sides)
  alpha5 = (new_center5 - sides)/(as.numeric(center_df[x,7]) - sides)
  alpha6 = (new_center6 - sides)/(as.numeric(center_df[x,7]) - sides)
  Rep4[[x]] <- round(new_center4 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha4))
  Rep5[[x]] <- round(new_center5 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha5))
  Rep6[[x]] <- round(new_center6 - (decay_matrix_subs*as.numeric(center_df[x,7])*alpha6))
}

#ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)
for(i in 1:length(de_half_cnv)){
  x <- de_half_cnv[i]
  new_center4 <- 1.5*as.numeric(center_df[x,4])
  new_center5 <- 1.5*as.numeric(center_df[x,5])
  new_center6 <- 1.5*as.numeric(center_df[x,6])
  Rep4[[x]] <- round(new_center4 - decay_matrix_subs*new_center4)
  Rep5[[x]] <- round(new_center5 - decay_matrix_subs*new_center5)
  Rep6[[x]] <- round(new_center6 - decay_matrix_subs*new_center6)
}

ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)

ALL_simulated <- array(unlist(ALL),dim=c(21,21,10000,6))
#ALL_simulated <- array(unlist(ALL),dim=c(31,31,10000,6))
ALL <- DelayedArray::DelayedArray(ALL_simulated)
Rep2[3]
spacings <- dim(ALL)
nBlocks <- 20
spacings[3] <- ceiling(dim(ALL)[3] / nBlocks)

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

grid <- RegularArrayGrid(dim(ALL), spacings)
grid
ans <-
  lapply(1:15, \(i) {
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

#Use direction function
##Function to get median of pixel parallel to diagnol at every MH distance
radial_par <- function(buffer, loop, inner){
  center <- buffer+1
  output <- c()
  output1 <- loop[center+inner,center+inner]
  output2 <- loop[center-inner,center-inner]
  output <- c(output1,output2)
  return(output)
}
dim(ALL)
#x <- decay_matrix_subs
#radial_par(buffer=15, decay_matrix_subs,inner=15)
ans <-
  lapply(1:10, \(i) {
    tmp <-
      blockApply(x = ALL,
                 FUN = \(x) apply(x,c(3,4),FUN = \(z) {
                   median(radial_par(buffer = buffer, loop=z, inner=i))
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
dim(ALL)
ans
observed <- ALL[11,11,,]
#observed <- as.matrix(ALL[16,16,,])
head(normalized)


## Isolate count matrix
loop <- data.frame(seq(1:10000))
loop1 <- cbind(loop,observed)
colnames(loop1) <- c("loop","WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
norm_MH = as.matrix(normalized) / exp(rowMeans(log(as.matrix(normalized))))
head(norm_MH)
to_remove <- which(!is.finite(rowSums(norm_MH)) == TRUE)
norm_MH <- norm_MH[-to_remove,]
observed <- observed[-to_remove,]
normalized <- normalized[-to_remove,]
new <- ALL[,,-to_remove,]
loops <- loop1[-to_remove,]
l <- subset(loops, loop %in% de)

dim(new)
dim(loops)
dim(observed)

##Add column names
colnames(observed) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
colnames(normalized) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
colnames(norm_MH) <- c("WT_1","WT_2","WT_3","TRT_1","TRT_2","TRT_3")
observed
head(normalized)
head(norm_MH)
##No normalization
colData <-
  do.call(rbind, strsplit(x = colnames(observed), split = "_")) |>
  as.data.frame(stringsAsFactors = TRUE) |>
  `colnames<-`(value = c("condition", "biorep"))
colData
observed <- abs(observed)
#observed1 <- as.matrix(observed[-1 ,])
dds <-
  DESeqDataSetFromMatrix(countData = round((observed)),
                         colData = colData,
                         design = ~ condition)
dds$condition <- factor(dds$condition, levels=c("WT","TRT"))
#default <- counts(dds,normalized = FALSE)
#default_sizefactor <- counts(dds, normalized =TRUE)
sizeFactors(dds)
#head(default_sizefactor)
## Run DEseq analysis

dds <- DESeq(dds)
res1 <-
  DESeq(dds) |>
  lfcShrink(coef = "condition_TRT_vs_WT", type="apeglm")
summary(res1)
file <- cbind(loops,res1)
head(file)
write.csv(file,"no_normalization_4-6.csv")
##With normalization
dds <-
  DESeqDataSetFromMatrix(countData = round(observed),
                         colData = colData,
                         design = ~ condition)

#normalizationFactors(dds) <- as.matrix(norm_MH)
#default <- counts(dds,normalized = FALSE)
#head(default)
#default_normfactor <- counts(dds, normalized =TRUE)
#head(default_normfactor)
#head(norm_MH)
dds <- DESeq(dds)
# perform enrichments
plotDispEsts(dds)
dds$condition <- factor(dds$condition, levels=c("WT","TRT"))
normalizationFactors(dds) <- as.matrix(norm_MH)

## Run DEseq analysis
res1 <-
#  DESeq(dds,fitType = 'mean')
  DESeq(dds) |>
  lfcShrink(coef = "condition_TRT_vs_WT")
summary(res1)
file1 <- cbind(loops,res1)
head(file1)
#write.csv(file2,"normalization_8-10.csv")
dim(subset(file, loop %in% de_half) |>
      subset(log2FoldChange >= 0.585) |>
      subset(padj < 0.1))
dim(subset(file, loop %in% de_half_cnv) |>
      subset(log2FoldChange >= 0.585) |>
      subset(padj < 0.1))
dim(subset(file1, loop %in% de_half) |>
    subset(log2FoldChange >= 0.585) |>
    subset(padj < 0.1))
dim(subset(file1, loop %in% de_half_cnv) |>
  subset(log2FoldChange >= 0.585) |>
  subset(padj < 0.1))
dim(subset(file, loop %in% de_half) |>
      subset(padj < 0.1))
dim(subset(file, loop %in% de_half_cnv) |>
      subset(padj < 0.1))
dim(subset(file1, loop %in% de_half) |>
      subset(padj < 0.1))
dim(subset(file1, loop %in% de_half_cnv) |>
      subset(padj < 0.1))

center_df[4155,]

dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst
res1 <-
  DESeq(dds,fitType='mean') |>
  lfcShrink(coef = "condition_TRT_vs_WT")
summary(res1)
