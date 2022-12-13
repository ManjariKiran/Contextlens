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
