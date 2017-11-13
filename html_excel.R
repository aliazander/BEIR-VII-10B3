#set this up for where you want to be working (working directory)
setwd("/Users/g-woloschak/Documents")

#run this if you've never run it before, but then you'll never need to again
install.packages('XML')
install.packages('plry')
install.packages('dplyr')

#run this every time you reopen R, but not every time you run a program
library(XML)
library(plyr)
library(dplyr)

#read in HTML and create table - copy and paste your url
df = readHTMLTable('file:///Users/g-woloschak/Downloads/homerResults.html', 
                   stringsAsFactors = FALSE)
#extract the first table
df <- df[[1]]

#change stupid names so they don't have special characters
names(df)[names(df)=="P-value"] <- "pvalue"
names(df)[names(df)=="log P-pvalue"] <- "lpvalue"
names(df)[names(df)=="Best Match/Details"] <- "details"

#select for the columns we are interested in
df <- select(df, pvalue, lpvalue, details)

#change p-value and log p-value from strings to numerics
df[,1] <- as.numeric(df[,1])
df[,2] <- as.numeric(df[,2])

#remove junk
df[-1] <- data.frame(lapply(df[-1], gsub, 
                 pattern = "More Information | Similar Motifs Found", 
                 replacement = "", 
                 fixed = TRUE))
#output file, name it whatever you want and it goes in your wd you set at the beginning
write.csv(df, file = "data_for_Dan.csv",
          row.names=TRUE)
