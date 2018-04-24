Factorial_reduce <- function(data){
  if (data==0){
  1
  } else
  reduce(c(1:data), function(x, y) x*y)
}
