#!/Users/vassiliadisdane/anaconda3/bin/R

options(warn=-1)
suppressMessages(library("tidyverse"))

args = commandArgs(TRUE)
path = file.path(args[1])
files <- list.files(path = path, pattern = "*_rawcounts.txt", recursive = T, include.dirs = FALSE, full.names = T)
message(paste("combining: ", files))
#print(files)
colnames <- c("Barcode")
for (i in files){
    names <- basename(i)
    #print(names)
    name <- strsplit(names, split = "_")[[1]][1]
    #print(name)
    colnames <- c(colnames, name)
}
#print(colnames)
input <- lapply(X = files, FUN = read.delim, header = F, sep = "\t")
input.merge <- input %>% reduce(left_join, by = "V1")
colnames(input.merge) <- colnames

write.table(input.merge, file = args[2], sep = "\t", quote = F, col.names=NA, row.names = T)
message("combine counts complete!")
