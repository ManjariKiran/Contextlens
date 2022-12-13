library(strawr)
library(tidyverse)
library(dbscan)
library(InteractionSet)
library(raster)
library(HiCcompare)
library(hictoolsr)


## Merge loops and convert to GInteractions

loops <- mergeBedpe(bedpeFiles = c("WT_5kbLoops.txt", "FS_5kbLoops.txt"), res = 10e3) |> as_ginteractions()
hicFiles <- c("GSM4259900_HEK_HiC_NUP_IDR_FS_A9_1_1_inter_30.hic")
l <- calcApa(bedpe = loops[1], hic = hicFiles[1], norm = "NONE", res = 10000, buffer = 10,filter = FALSE)
l
mh_angle <- function(buffer, loop, direction){
  m=(buffer*2)+1
  center <- buffer+1
  center1 <- (m+1)/2
  M <- matrix(data=NA,nrow=m,ncol=m)
  for(i in 1:center1){
    for(j in 1:center1){
      base = center-i
      height = center-j
      angle = (round(atan(base/height)*(180/pi)))
      M[i,j]<-angle
    }
  }
  for(i in 1:m){
    for(j in center:m){
      base = i-center
      height = j-center
      angle = 180+(round(atan(base/height)*(180/pi)))
      M[i,j]<-angle
    }
  }
  for(i in center:m){
    for(j in 1:buffer){
      base = i-center
      height = j-center
      angle = 360+(round(atan(base/height)*(180/pi)))
      M[i,j]<-angle
    }
  }
  val <- M %in% direction
  v <- which(val == "TRUE")
  new <- loop[v]
  return(new)
}

mh_angle(3,l,direction=90:180)


l <- calcApa(bedpe = loops[1], hic = hicFiles[1], norm = "NONE", res = 10000, buffer = 4,filter = FALSE)

donut <- function(buffer, loop,inner){
  m=(buffer*2)+1
  center <- buffer+1
  up=center-inner
  down=center+inner
  M <- matrix(data=1,nrow=m,ncol=m)
  for(i in 1:m){
    M[i,1]<-0
    M[i,m]<-0
    M[i,center]<-0
  }
  for(j in 1:m){
    M[1,j]<-0
    M[m,j]<-0
    M[center,j]<-0
  }
  for(i in up:down){
    for(j in up:down){
      M[i,j]<-0
    }
  }
  inner_val <- which(M == 1)
  new <- loop[inner_val]
  return(new)
}

donut(4,l,1)
l
