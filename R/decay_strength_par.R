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

radial_par <- function(buffer, loop, inner){
  center <- buffer+1
  output <- c()
  output1 <- loop[center+inner,center+inner]
  output2 <- loop[center-inner,center-inner]
  output <- c(output1,output2)
  return(output)
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

val1 <- c()
foo1 <- data.frame()
output1 <- data.frame()
for(i in 1:dim(countMatrix_obs1)[3]){
  val1 <- append(val1,countMatrix_obs1[11,11,i])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs1[,,i])
  foo1[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output1 <- cbind(val1,foo1)
}

val2 <- c()
foo2 <- data.frame()
output2 <- data.frame()
for(i in 1:dim(countMatrix_obs2)[3]){
  val2 <- append(val2,countMatrix_obs2[11,11,i])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs2[,,i])
    foo2[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output2 <- cbind(val2,foo2)
}

val3 <- c()
foo3 <- data.frame()
output3 <- data.frame()
for(i in 1:dim(countMatrix_obs3)[3]){
  val3 <- append(val3,countMatrix_obs3[11,11,i])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs3[,,i])
    foo3[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output3 <- cbind(val3,foo3)
}

val4 <- c()
foo4 <- data.frame()
output4 <- data.frame()
for(i in 1:dim(countMatrix_obs4)[3]){
  val4 <- append(val4,countMatrix_obs4[11,11,i])
  for(j in 1:10){
    loop=as.matrix(countMatrix_obs4[,,i])
    foo4[i,j] <- median(radial_par(buffer = buffer, loop=loop, inner=j))
  }
  output4 <- cbind(val4,foo4)
}


#Convert data by setting 100 kb as 0 (subtract 100kb value from all values) and center pixel as 1 (divide all values to center pixel)
zero1=data.frame()
zero2=data.frame()
zero3=data.frame()
zero4=data.frame()
dim(countMatrix_obs1)
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
data
