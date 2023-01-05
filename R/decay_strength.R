#Program to show that decay remain same for different Strength, distance and loop strength
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

## Select different Strength
loops <- binPairs(x=loops, binSize=10e03)

head(loops)
## Modify seqlevels of loops to match hicFiles
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'


## Expand pixels to matrices and extract
loopCounts <-
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(binSize=10e03,
                  files=hicFiles[5],
                  norm=norm,
                  matrix=matrix)

loops
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


countMatrix_obs <-
  assay(loopCounts) |>
  aperm(c(3,4,1,2))

x <- (countMatrix_obs[11,11,,])
quantile(x)
as.numeric(quantile(x,0.25))
as.numeric(quantile(x,0.50))
as.numeric(quantile(x,0.75))
a <- which(x <= as.numeric(quantile(x,0.25)))
b <- which(x <= as.numeric(quantile(x,0.50)) & (x > as.numeric(quantile(x,0.25))))
c <- which(x <= as.numeric(quantile(x,0.75)) & (x > as.numeric(quantile(x,0.50))))
d <- which(x > as.numeric(quantile(x,0.75)))
countMatrix_obs1 <- (countMatrix_obs[,,a,])
countMatrix_obs2 <- (countMatrix_obs[,,b,])
countMatrix_obs3 <- (countMatrix_obs[,,c,])
countMatrix_obs4 <- (countMatrix_obs[,,d,])
dim(countMatrix_obs1)
dim(countMatrix_obs2)
dim(countMatrix_obs3)
dim(countMatrix_obs4)
spacings <- dim(countMatrix_obs1)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs1)[1] / nBlocks)
grid <- RegularArrayGrid(dim(countMatrix_obs1), spacings)
ans_l1 <-
  lapply(0:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs1,
                 FUN = \(x) apply(x,c(3),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })
spacings <- dim(countMatrix_obs2)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs2)[1] / nBlocks)
grid <- RegularArrayGrid(dim(countMatrix_obs2), spacings)
ans_l2 <-
  lapply(0:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs2,
                 FUN = \(x) apply(x,c(3),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })
spacings <- dim(countMatrix_obs3)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs3)[1] / nBlocks)
grid <- RegularArrayGrid(dim(countMatrix_obs3), spacings)
ans_l3 <-
  lapply(0:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs3,
                 FUN = \(x) apply(x,c(3),FUN = \(z) {
                   median(mh_index(buffer = buffer, loop=z, inner=i, exclude= 0))
                 }),
                 grid = grid,
                 verbose = TRUE,
                 BPPARAM = MulticoreParam())
    do.call(rbind, tmp)
  })
spacings <- dim(countMatrix_obs4)
nBlocks <- 20
spacings[3] <- ceiling(dim(countMatrix_obs4)[1] / nBlocks)
grid <- RegularArrayGrid(dim(countMatrix_obs4), spacings)
ans_l4 <-
  lapply(0:10, \(i) {
    tmp <-
      blockApply(x = countMatrix_obs4,
                 FUN = \(x) apply(x,c(3),FUN = \(z) {
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
dim(countMatrix_obs1)
for(j in 1:11){
  for(i in 1:dim(countMatrix_obs1)[3]){
    a <- ans_l1[[1]][[i]]-ans_l1[[11]][[i]]
    b <- ans_l1[[j]][[i]]-ans_l1[[11]][[i]]
    zero1[i,j] <- round((b/a),2)
  }
}

for(j in 1:11){
  for(i in 1:dim(countMatrix_obs2)[3]){
    a <- ans_l2[[1]][[i]]-ans_l2[[11]][[i]]
    b <- ans_l2[[j]][[i]]-ans_l2[[11]][[i]]
    zero2[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
  for(i in 1:dim(countMatrix_obs3)[3]){
    a <- ans_l3[[1]][[i]]-ans_l3[[11]][[i]]
    b <- ans_l3[[j]][[i]]-ans_l3[[11]][[i]]
    zero3[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
  for(i in 1:dim(countMatrix_obs4)[3]){
    a <- ans_l4[[1]][[i]]-ans_l4[[11]][[i]]
    b <- ans_l4[[j]][[i]]-ans_l4[[11]][[i]]
    zero4[i,j] <- round((b/a),2)
  }
}

colnames(zero1) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero2) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero3) <- c("0","1","2","3","4","5","6","7","8","9","10")
colnames(zero4) <- c("0","1","2","3","4","5","6","7","8","9","10")

Group <- factor(0:10)
normalized_1 <- na.omit(zero1)
normalized_2 <- na.omit(zero2)
normalized_3 <- na.omit(zero3)
normalized_4 <- na.omit(zero4)

Mean1<- as.vector(colMeans(normalized_1,na.rm=TRUE))
Mean2<- as.vector(colMeans(normalized_2,na.rm=TRUE))
Mean3<- as.vector(colMeans(normalized_3,na.rm=TRUE))
Mean4<- as.vector(colMeans(normalized_4,na.rm=TRUE))

Res1 <- c("First")
Res2 <- c("Second")
Res3 <- c("Third")
Res4 <- c("Fourth")

###Plotting ########
df1 <- data.frame(Group,Mean1,Res1)
colnames(df1) <- c("Group","Mean","Strength")
df2 <- data.frame(Group,Mean2,Res2)
colnames(df2) <- c("Group","Mean","Strength")
df3 <- data.frame(Group,Mean3,Res3)
colnames(df3) <- c("Group","Mean","Strength")
df4 <- data.frame(Group,Mean4,Res4)
colnames(df4) <- c("Group","Mean","Strength")
data <- rbind(df1,df2,df3,df4)
data$Strength
library(ggplot2)
library(forcats)
colPalette <- c("#8856a7","#8c96c6","#9ebcda","#bfd3e6")
colPalette <- c("#bfd3e6","#9ebcda","#8c96c6","#8856a7")
p <- ggplot(data,aes(Group,Mean,group=fct_inorder(Strength)))+geom_line(aes(color=fct_inorder(Strength)))+
  geom_point(aes(color=Strength))+xlab("Distance from the center") + ylab("Signal relative to the center")+theme_classic() +scale_y_continuous(limits = c(0, 1)) + scale_color_manual(values=colPalette)
p
p + theme(
  axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),
  axis.title.x = element_text(color="grey", size=20, face="bold"),
  axis.title.y = element_text(color="grey", size=20, face="bold"),
  legend.position="right",legend.text = element_text(size=16),legend.title = element_text(size=16)
)

