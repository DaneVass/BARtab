#!/usr/bin/env Rscript

options(warn=-1)
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("purrr"))

args = commandArgs(TRUE)
files = args[1:length(args)-1]
outname = args[length(args)]

print(args)
print(files)
print(outname)

#path = file.path(args[1])
#files <- list.files(path = path, pattern = "*_rawcounts.txt", recursive = T, include.dirs = FALSE, full.names = T)
files <- files[grep("*_rawcounts.txt", files)]
print(files)

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

write.table(input.merge, file = outname, sep = "\t", quote = F, col.names=NA, row.names = T)
message("combine counts complete!")
