## Part 1: Factorial Function ----

## Factorial_loop: a version that computes the factorial of an integer using looping (such as a for loop) ----
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

## Factorial_reduce: a version that computes the factorial using the reduce() function in the purrr package ----
Factorial_reduce <- function(n){
  if (n==0) 1 else purrr::reduce(c(1:n), function(x, y) x*y)
}

## Factorial_func: a version that uses recursion to compute the factorial ----                                  
Factorial_func <- function(n) {
  if (n == 0) 1 else n * Factorial_func(n-1)
} 
                                 
## Factorial_mem: a version that uses memoization to compute the factorial ----
                                 
mem <- function(){
  
  out <- 1
  
  Factorial_mem <- function(x){
    
    if(x < 0) stop()
    
    if (x == 0 | x == 1) return(1)
    
    if (length(out) < x) out <<- `length<-`(out, x)
    
    out[x] <<- x * factorial(x-1)
    out[x]
  }
  Factorial_mem
}

Factorial_mem <- mem()
