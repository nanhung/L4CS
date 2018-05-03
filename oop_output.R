df<-read.csv("MIE.csv")
class(df)
names(df)

make_LD <- function(df){
  if (class(df)=="data.frame") x <- structure(df, class = "LongitudinalData")
  else stop ("The input should be a data frame!")
  str(x)
  invisible(x)
}


subject <- function(x, id) UseMethod("subject")
visit <- function(x, id) UseMethod("visit")
room <- function(x, id) UseMethod("room")

subject.LongitudinalData <- function(x, i) {
  
  index <- which(x$id %in% i)
  x <- lapply(x, function(x) x[index])
  
  structure(x, class = "LongitudinalData")
}

############
library(pryr)

x<-make_LD(df)


ftype(make_LD)
ftype(subject)
