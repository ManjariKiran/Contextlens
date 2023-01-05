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

countMatrix_obs1 <-
  assay(loopCounts1) |>
  aperm(c(3,4,1,2))
dim(countMatrix_obs1)
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

radial_par <- function(buffer, loop, inner){
  center <- buffer+1
  output <- c()
  output1 <- loop[center+inner,center+inner]
  output2 <- loop[center-inner,center-inner]
  output <- c(output1,output2)
  return(output)
}

val1 <- c()
foo1 <- data.frame()
output1 <- data.frame()
for(i in 1:dim(countMatrix_obs1)[3]){
  val1 <- append(val1,countMatrix_obs1[11,11,i,1])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs1[,,i,1])
    foo1[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output1 <- cbind(val1,foo1)
}

val2 <- c()
foo2 <- data.frame()
output2 <- data.frame()
for(i in 1:dim(countMatrix_obs2)[3]){
  val2 <- append(val2,countMatrix_obs2[11,11,i,1])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs2[,,i,1])
    foo2[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output2 <- cbind(val2,foo2)
}

val3 <- c()
foo3 <- data.frame()
output3 <- data.frame()
for(i in 1:dim(countMatrix_obs3)[3]){
  val3 <- append(val3,countMatrix_obs3[11,11,i,1])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs3[,,i,1])
    foo3[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output3 <- cbind(val3,foo3)
}

val4 <- c()
foo4 <- data.frame()
output4 <- data.frame()
for(i in 1:dim(countMatrix_obs4)[3]){
  val4 <- append(val4,countMatrix_obs4[11,11,i,1])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs4[,,i,1])
    foo4[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output4 <- cbind(val4,foo4)
}

val5 <- c()
foo5 <- data.frame()
output5 <- data.frame()
for(i in 1:dim(countMatrix_obs5)[3]){
  val5 <- append(val5,countMatrix_obs5[11,11,i,1])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs5[,,i,1])
    foo5[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output5 <- cbind(val5,foo5)
}


#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero1=data.frame()
zero2=data.frame()
zero3=data.frame()
zero4=data.frame()
zero5=data.frame()


for(j in 1:11){
  for(i in 1:dim(countMatrix_obs1)[3]){
    a <- output1[i,1]-output1[i,11]
    b <- output1[i,j]-output1[i,11]
    zero1[i,j] <- round((b/a),2)
  }
}

for(j in 1:11){
  for(i in 1:dim(countMatrix_obs2)[3]){
    a <- output2[i,1]-output2[i,11]
    b <- output2[i,j]-output2[i,11]
    zero2[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
  for(i in 1:dim(countMatrix_obs3)[3]){
    a <- output3[i,1]-output3[i,11]
    b <- output3[i,j]-output3[i,11]
    zero3[i,j] <- round((b/a),2)
  }
}
for(j in 1:11){
  for(i in 1:dim(countMatrix_obs4)[3]){
    a <- output4[i,1]-output4[i,11]
    b <- output4[i,j]-output4[i,11]
    zero4[i,j] <- round((b/a),2)
  }
}

for(j in 1:11){
  for(i in 1:dim(countMatrix_obs5)[3]){
    a <- output5[i,1]-output5[i,11]
    b <- output5[i,j]-output5[i,11]
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
colPalette <- c("#bfd3e6","#9ebcda","#8c96c6","#8856a7","#810f7c")
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
