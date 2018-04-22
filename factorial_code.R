Factorial_loop <- function(x){
  if (x==0){
    y <- 1
  } else {
    y<-1
    for (i in 1:x){
      y <-y*((1:x)[i])
    }
  }
  print(y)
}
