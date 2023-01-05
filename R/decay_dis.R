#Program to show that decay remain same for different resolution, distance and loop strength
library(strawr)
library(tidyverse)
library(dbscan)
library(InteractionSet)
library(raster)
library(HiCcompare)
library(hictoolsr)

## List loop files
loopFiles <- list.files(path = "~/contextlens/inst/extdata",pattern=".txt$", full.names = TRUE)
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
head(loops)
## Modify seqlevels of loops to match hicFiles
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

head(loops)
a <- pairdist(loops) > 0
l1 <- loops[which(a[a == TRUE])]
a <- pairdist(loops) > 50000
l2 <- loops[which(a[a == TRUE])]
a <- pairdist(loops) > 100000
l3 <- loops[which(a[a == TRUE])]
a <- pairdist(loops) > 200000
l4 <- loops[which(a[a == TRUE])]
a <- pairdist(loops) > 300000
l5 <- loops[which(a[a == TRUE])]

length(l1)
length(l2)
length(l3)
length(l4)
length(l5)

## Expand pixels to matrices and extract
loopCounts1 <-
  pixelsToMatrices(x=l1, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)
loopCounts2 <-
  pixelsToMatrices(x=l2, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)
loopCounts3 <-
  pixelsToMatrices(x=l3, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)
loopCounts4 <-
  pixelsToMatrices(x=l4, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)

loopCounts5 <-
  pixelsToMatrices(x=l5, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)

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


countMatrix_obs1 <-
  assay(loopCounts1) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs)
countMatrix_obs2 <-
  assay(loopCounts2) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs2)
countMatrix_obs3 <-
  assay(loopCounts3) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs3)
countMatrix_obs4 <-
  assay(loopCounts4) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs4)
countMatrix_obs5 <-
  assay(loopCounts5) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs5)

spacings <- dim(countMatrix_obs)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs)[3] / nBlocks)
grid <- RegularArrayGrid(dim(countMatrix_obs1), spacings)

ans_l1 <-
  lapply(0:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs1,
                 FUN = \(x) apply(x,c(3,4),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })

#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero1=data.frame()
zero2=data.frame()
zero3=data.frame()
zero4=data.frame()
zero5=data.frame()

for(j in 1:11){
for(i in 1:length(l1)){
    a <- ans_l1[[1]][[i]]-ans_l1[[11]][[i]]
    b <- ans_l1[[j]][[i]]-ans_l1[[11]][[i]]
    zero1[i,j] <- round((b/a),2)
  }
}

for(i in 1:length(loops)){
  for(j in 1:11){
    a <- normalized[i,1]-normalized[i,11]
    b <- normalized[i,j]-normalized[i,11]
    zero[i,j] <- round(((b)/(a)),2)
  }
}

for(j in 1:11){
    for(i in 1:length(l2)){
    a <- ans_l2[[1]][[i]]-ans_l2[[11]][[i]]
    b <- ans_l2[[j]][[i]]-ans_l2[[11]][[i]]
    zero2[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
for(i in 1:length(l3)){
    a <- ans_l3[[1]][[i]]-ans_l3[[11]][[i]]
    b <- ans_l3[[j]][[i]]-ans_l3[[11]][[i]]
    zero3[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
for(i in 1:length(l4)){
    a <- ans_l4[[1]][[i]]-ans_l4[[11]][[i]]
    b <- ans_l4[[j]][[i]]-ans_l4[[11]][[i]]
    zero4[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
for(i in 1:length(l5)){
    a <- ans_l5[[1]][[i]]-ans_l5[[11]][[i]]
    b <- ans_l5[[j]][[i]]-ans_l5[[11]][[i]]
    zero5[i,j] <- round((b/a),2)
  }
}

colnames(zero1) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero2) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero3) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero4) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero5) <- c("0","1","2","3","4","5","6","7","8","9","10")

Group <- factor(0:10)
normalized_1 <- na.omit(zero1)
normalized_2 <- na.omit(zero2)
normalized_3 <- na.omit(zero3)
normalized_4 <- na.omit(zero4)
normalized_5 <- na.omit(zero5)
Mean1<- as.vector(colMeans(normalized_1,na.rm=TRUE))
Mean2<- as.vector(colMeans(normalized_2,na.rm=TRUE))
Mean3<- as.vector(colMeans(normalized_3,na.rm=TRUE))
Mean4<- as.vector(colMeans(normalized_4,na.rm=TRUE))
Mean5<- as.vector(colMeans(normalized_5,na.rm=TRUE))
normalized_1
Res1 <- c("0")
Res2 <- c("50")
Res3 <- c("100")
Res4 <- c("200")
Res5 <- c("300")
Mean1
###Plotting ########
df1 <- data.frame(Group,Mean1,Res1)
colnames(df1) <- c("Group","Mean","Distance")
df2 <- data.frame(Group,Mean2,Res2)
colnames(df2) <- c("Group","Mean","Distance")
df3 <- data.frame(Group,Mean3,Res3)
colnames(df3) <- c("Group","Mean","Distance")
df4 <- data.frame(Group,Mean4,Res4)
colnames(df4) <- c("Group","Mean","Distance")
df5 <- data.frame(Group,Mean5,Res5)
colnames(df5) <- c("Group","Mean","Distance")

data <- rbind(df1,df2,df3,df4,df5)
data
library(ggplot2)
colPalette <- c("#810f7c","#8856a7","#8c96c6","#9ebcda","#bfd3e6")
p <- ggplot(data,aes(Group,Mean,group=fct_inorder(Distance)))+geom_line(aes(color=fct_inorder(Distance)))+
  geom_point(aes(color=Distance))+xlab("Distance from the center") + ylab("Signal relative to the center")+theme_classic() +scale_y_continuous(limits = c(0, 1)) + scale_color_manual(values=colPalette)
p
p + theme(
  axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),
  axis.title.x = element_text(color="grey", size=20, face="bold"),
  axis.title.y = element_text(color="grey", size=20, face="bold"),
  legend.position="right",legend.text = element_text(size=16),legend.title = element_text(size=16)
)
normalized_1
