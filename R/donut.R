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

#donut(4,l,1)
