library(strawr)
library(devtools)
library(mariner)
library(SummarizedExperiment)
library(DelayedArray)
library(BiocParallel)
library(InteractionSet)

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

## Expand pixels to matrices and extract
loopCounts <-
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(binSize=res,
                  files=hicFiles[1],
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

observed <- as.matrix(countMatrix_obs[11,11,,])
head(observed)

decay <- c()
x <- c()
y <- c()
z <- c()
for(i in length(loops)){  #for all the loops in the example file
  center <- countMatrix_obs[11,11,i,1]  #Extract the center pixel
  #diff <- (center - ans[[1]])/center  #difference between center and MH distance 1 pixel
  #x <- as.vector(colMeans(diff)) #Mean of the difference
  for(j in 1:10){   #Now for remining MH distances
    #  start <- ans[[j]]  #Start with MH distance 1
    #  k <- j+1
    #  end <- ans[[k]]   #End with the next MH distance
    diff1 <- (center - ans[[j]])/center   #Calculate the difference between start and end
    y <-colMeans(diff1)    #Mean of the difference at each step
    decay <- as.vector(append(decay,y)) #Append the vector
    #  decay <- c(x,z)
  }
}
decay

mh_index <- function(buffer, loop_signal, decay){
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
  M[center,center] <- loop_signal
  for(i in 1:8){
    MH_val <- which(M == i)
    val <- loop_signal*decay[i]
    M[MH_val] <- val
  }
  return(M)
}

