radial_direction <- function(buffer, loop, direction, inner){
  center <- buffer+1
  if(direction == 3){
    output <- loop[center+inner,center+inner]
  }
  if(direction == 1){
    output <- loop[center-inner,center-inner]
  }
  if(direction == 2){
    output <- loop[center-inner,center+inner]
  }
  if(direction == 4){
    output <- loop[center+inner,center-inner]
  }

  return(output)
}

radial_par <- function(buffer, loop, inner){
  center <- buffer+1
  output <- c()
    output1 <- loop[center+inner,center+inner]
    output2 <- loop[center-inner,center-inner]
    output <- c(output1,output2)
  return(output)

}
l1 <- as.matrix(countMatrix_obs1[,,1])
l1
radial_par(10,l1,1)
l1
