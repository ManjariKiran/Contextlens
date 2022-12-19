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
original <- 100 - (100*decay_matrix_subs)
sides <- median((c(original[1,],original[,1],original[m,],original[,m])),na.rm=TRUE)
alpha = (200 - sides)/(100 - sides)

decay_matrix_subs
TRT <- round(200 - (100*decay_matrix_subs*alpha))
WT <- 100 - (100*decay_matrix_subs)

WT <- decay_matrix*100
TRT <- decay_matrix*200
normalized1=c()
normalized2=c()
for(i in 0:10){
  new1 <- median(mh_index(buffer = 10, loop = WT , inner = i,exclude=0))
  new2 <- median(mh_index(buffer = 10, loop = TRT , inner = i,exclude=0))
  normalized1 <- append(normalized1,new1)
  normalized2 <- append(normalized2,new2)
}


observed <- as.matrix(countMatrix_obs[11,11,,])
as.vector(countMatrix_obs[11,11,1,])

##Calculate mean and sd of center pixel
observed <- as.matrix(countMatrix_obs[11,11,,])
colnames(observed) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2")
dim(observed)

loop_mean <- rowMeans(observed)
loop_sd <- rowSds(observed)
observed_stat <- cbind(observed,round(loop_mean))
colnames(observed_stat) <- c("FS_1_1","FS_1_2","FS_2_1","FS_2_2","mean")
number <- sort(unique(observed_stat[,5]))
ref <- data.frame("mean","sd")
colnames(ref) <- c("mean","sd")
for(i in 1:length(number)){
  x <- mean(observed_stat[which(observed_stat[,5] == number[i]),])
  y <- sd(observed_stat[which(observed_stat[,5] == number[i]),])
  ref[i,] <- cbind(number[i],round(y,3))
}
head(observed_stat)

number
##Print random numbers for given mean and std
center_df <- data.frame("rep1","rep2","rep3","rep4","rep5","rep6")
colnames(center_df) <- c("rep1","rep2","rep3","rep4","rep5","rep6")
for(i in 1:10000){
   random <-  sample(number, 1, replace = TRUE)
   n_sd <- as.numeric(ref[which(ref$mean == random),]$sd)
   center_df[i,] <- round(rnorm(n=6,mean = random, sd= n_sd),0)
}

##Set these as center pixel and create local background
head(center_df)
center_df2 <- center_df[,4:6]
head(center_df2)
##Create 10% fold change loops
de <- sample(1:10000,size = 1000, replace=F)
length(de)
write.csv(de,"total_de")
de_half <- sample(de,size=400,replace=F)
length(de_half)
write.csv(de_half,"total_de_half")
de_rem <- de[!de %in% de_half]
length(de_rem)
de_two <- sample(de_rem,size=300,replace=F)
length(de_two)
write.csv(de_two,"total_de_two")
de_rem <- de_rem[!de_rem %in% de_two]
length(de_rem)
de_three <- sample(de_rem,size=200,replace=F)
length(de_three)
write.csv(de_three,"total_de_three")
de_rem <- de_rem[!de_rem %in% de_three]
length(de_rem)
de_four <- de_rem
write.csv(de_four,"total_de_four")
de_four
head(center_df2)
#de_four <-(center_df2[de_four,][1:3])
#de_four[1:3] <- apply(de_four[1:3],2,
#                    function(x) as.numeric(as.character(x)))
#de_three <-(center_df2[de_three,][1:3])
#de_three[1:3] <- apply(de_three[1:3],2,
#                      function(x) as.numeric(as.character(x)))
#de_two <-(center_df2[de_two,][1:3])
#de_two[1:3] <- apply(de_two[1:3],2,
#                      function(x) as.numeric(as.character(x)))
#de_half <-(center_df2[de_half,][1:3])
#de_half[1:3] <- apply(de_half[1:3],2,
#                     function(x) as.numeric(as.character(x)))

#de_half
#Remaining <-center_df2[-de,][1:3]
#center_df2 <- rbind(Remaining,de_half,de_two,de_three,de_four)
#treatment <- as.data.frame(center_df2)
#dim(treatment)
#center_df2 <- treatment[order(as.numeric(row.names(center_df2))),]
#head(center_df)
#center_df <- center_df[,1:3]

Rep1 <- vector("list",length = 10000)
Rep2 <- vector("list",length = 10000)
Rep3 <- vector("list",length = 10000)
Rep4 <- vector("list",length = 10000)
Rep5 <- vector("list",length = 10000)
Rep6 <- vector("list",length = 10000)
m=(buffer*2)+1
for(i in 1:10000){
  M1 <- matrix(data=NA,nrow=m,ncol=m)
  M2 <- matrix(data=NA,nrow=m,ncol=m)
  M3 <- matrix(data=NA,nrow=m,ncol=m)
  M4 <- matrix(data=NA,nrow=m,ncol=m)
  M5 <- matrix(data=NA,nrow=m,ncol=m)
  M6 <- matrix(data=NA,nrow=m,ncol=m)
      M1 <- round(as.numeric(center_df[i,1]) - decay_matrix_subs*as.numeric(center_df[i,1]))
      M2 <- round(as.numeric(center_df[i,2]) - decay_matrix_subs*as.numeric(center_df[i,2]))
      M3 <- round(as.numeric(center_df[i,3]) - decay_matrix_subs*as.numeric(center_df[i,3]))
      M4 <- round(as.numeric(center_df[i,4]) - decay_matrix_subs*as.numeric(center_df[i,4]))
      M5 <- round(as.numeric(center_df[i,5]) - decay_matrix_subs*as.numeric(center_df[i,5]))
      M6 <- round(as.numeric(center_df[i,6]) - decay_matrix_subs*as.numeric(center_df[i,6]))
Rep1[[i]] <- M1
Rep2[[i]] <- M2
Rep3[[i]] <- M3
Rep4[[i]] <- M4
Rep5[[i]] <- M5
Rep6[[i]] <- M6
}
ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)
de_two
Rep1[[5913]]
for(i in 1:length(de_two)){
  x <- de_two[i]
  sides <- median((c(Rep1[[x]][1,],Rep1[[x]][,1],Rep1[[x]][m,],Rep1[[x]][,m])),na.rm=TRUE)
  new_center4 <- 2*as.numeric(center_df[x,4])
  new_center5 <- 2*as.numeric(center_df[x,5])
  new_center6 <- 2*as.numeric(center_df[x,6])
  alpha4 = (new_center4 - sides)/(as.numeric(center_df[x,1]) - sides)
  alpha5 = (new_center5 - sides)/(as.numeric(center_df[x,2]) - sides)
  alpha6 = (new_center6 - sides)/(as.numeric(center_df[x,3]) - sides)
  Rep4[[x]] <- round(new_center4 - (decay_matrix_subs*as.numeric(center_df[x,1])*alpha4))
  Rep5[[x]] <- round(new_center5 - (decay_matrix_subs*as.numeric(center_df[x,2])*alpha5))
  Rep6[[x]] <- round(new_center6 - (decay_matrix_subs*as.numeric(center_df[x,3])*alpha6))
}
Rep4[[5913]]
ALL <- list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6)
ALL_simulated <- array(unlist(ALL),dim=c(21,21,10000,6))
ALL <- DelayedArray::DelayedArray(ALL_simulated)
dim(ALL)
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
l <- subset(loops, loop %in% de_two)

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

# perform enrichments

dds$condition <- factor(dds$condition, levels=c("WT","TRT"))
normalizationFactors(dds) <- as.matrix(norm_MH)

resultsNames(dds)
## Run DEseq analysis
res1 <-
  DESeq(dds,fitType='local') |>
  lfcShrink(coef = "condition_TRT_vs_WT", type="apeglm")
summary(res1)
file1 <- cbind(loops,res1)
head(file1)
write.csv(file1,"normalization.csv")
subset(file1, loop %in% de_two) |>
    subset(log2FoldChange > 1) |>
    subset(padj < 0.1)




