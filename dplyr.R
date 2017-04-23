# 0423 Hadley Wickham's "dplyr" tutorial at useR 2014 (1/2)
# https://www.youtube.com/watch?v=8SGif63VW6E&t=1593s
rm(list = ls()) 

if(!require(dplyr)) {
  install.packages("dplyr"); require(dplyr)
}

if(!require(hflights)) {install.packages("hflights"); require(hflights)}

data(hflights)
head(hflights)
View(hflights)
ls(hflights)

flights <- tbl_df(hflights)

hflights
flights

 
class(flights)

sfo <- filter(flights, Dest == "SFO")
sfo

sfooak <- filter(flights, Dest %in% c("SFO","OAK"))

filter(flights, Month == 1)
filter(flights, DepDelay > 60)
filter(flights, ArrDelay > 2 * DepDelay)

select(flights, ArrDelay:DepDelay) 
select(flights, ArrDelay:DepDelay) 
select(flights, contains("Delay"))
select(flights, ends_with("Delay"))
arrange(flights, Year, Month, DayofMonth)
MaxDepDelay <- arrange(flights, desc(DepDelay))
flights<-mutate(flights, speed = Distance / (AirTime/60))
MaxSpeed <- arrange(flights, desc(speed))
