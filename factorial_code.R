Factorial_loop <- function(n){
  if (n==0){
    y <- 1
  } else {
    y <- 1
    for (i in 1:n){
      y <- y * ((1:n)[i])
    }
  }
  print(y)
}

Factorial_reduce <- function(n){
  if (n==0) 1 else purrr::reduce(c(1:n), function(x, y) x*y)
}

Factorial_func <- function(n) {
  if (n == 0) 1 else n * Factorial_func(n-1)
} 
