df<-read.csv("MIE.csv")
class(df)
names(df)

make_LD <- function(df){
  if (class(df)=="data.frame") x <- structure(df, class = "LongitudinalData")
  else stop ("The input have to be a data frame")
}

subject <- function(x, i) UseMethod("subject")
visit <- function(x, i) UseMethod("visit")
room <- function(x, i) UseMethod("room")

subject.LongitudinalData <- function(x, i) {

 object <- deparse(substitute(x))
  
  if (all(i %in% x$id)) {
    subject <- which(x$id %in% i)
    for(j in 1:length(x)){
      x[[j]] <- x[[j]][subject]
    }
  } else stop ("Not avaliable id")

  assign(object, x, parent.frame())  
}

visit.LongitudinalData <- function(x, i) {
  
  object <- deparse(substitute(x))
  
  if (all(i %in% x$visit)) {
    visit <- which(x$visit %in% i)
    for(j in 1:length(x)){
      x[[j]] <- x[[j]][visit]
    }
  } else stop ("Not avaliable number")
  
  assign(object, x, parent.frame())  
  
}

room.LongitudinalData <- function(x, i) {
  
  object <- deparse(substitute(x))
  
  if (all(i %in% x$room)) {
    room <- which(x$room %in% i)
    for(j in 1:length(x)){
      x[[j]] <- x[[j]][room]
    }
  } else stop ("Not avaliable room")
  
  assign(object, x, parent.frame())  
  
}

print.LongitudinalData <- function(x) {
  cat("\nLongitudinalData: n =", length(x$id), "\n")
  cat("\nID [", class(x$id),  "]\n")
  print(table(x$id))
  cat("\n")
  cat("\nVisit [", class(x$visit), "]\n")
  print(table(x$visit))
  cat("\n")
  cat("\nRoom [", class(x$room), "]\n")
  print(table(x$room))
  cat("\n")
  cat("\nValue [", class(x$value), "]\n")
  print(summary(x$value))
  cat("\n")
  cat("\nTimepoint [", class(x$timepoint), "]\n")
  print(summary(x$timepoint))
}


summary.LongitudinalData <- function(x) {
  cat("\nLongitudinalData: n =", length(x$id), "\n")
  cat("\nID\n")
  cat(unique(x$id), sep = ", ")
  cat("\n")
  cat("\nVisit\n")
  cat(unique(x$visit), sep = ", ")
  cat("\n")
  cat("\nRoom\n")
  cat( paste(unique(x$room)), sep = ", ")
  cat("\n")
  cat("\nValue\n")
  cat(min(x$value), "-",max(x$value), "(min-max)")
  cat("\n")
  cat("\nTimepoint\n")
  cat(min(x$timepoint), "-",max(x$timepoint), "(min-max)\n")
}



############
library(pryr)
library(dplyr)
ftype(make_LD)
ftype(subject)


x<-make_LD(df) 

make_LD(df) %>% summary

make_LD(df) %>% print # All data

make_LD(df) %>% subject(14) %>% print # subject 14

make_LD(df) %>% subject(14) %>% summary # subject 14

make_LD(df) %>% visit(1) %>% print # visit 1

make_LD(df) %>% room("den") %>% print # room den

make_LD(df) %>% subject(c(14,20)) %>% visit(c(1,2)) %>% print 

make_LD(df) %>% subject(c(14,20)) %>% visit(c(1,2)) %>% room(c("den", "bedroom")) %>% print




