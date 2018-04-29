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

Factorial_reduce <- function(data){
  library(purrr)
  if (data==0){
  1
  } else
  reduce(c(1:data), function(x, y) x*y)
}
         
Factorial_func <- function(x) {
  if (x == 0) {
    1
  }
  else {
    x * Factorial_func(x-1)
  }
}         
