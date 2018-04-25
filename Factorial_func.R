Factorial_func <- function(x) {
  if (x == 0) {
    1
  }
  else {
    x * Factorial_func(x-1)
  }
}
