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

#mh_angle(3,l,direction=90:180)
