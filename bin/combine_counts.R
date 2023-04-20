#!/usr/bin/env Rscript

options(warn = -1)
suppressMessages(library("dplyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("purrr"))

args <- commandArgs(TRUE)
files <- args[1:length(args) - 1]
outname <- args[length(args)]

print(files)
print(outname)

colnames <- c("Barcode")
for (i in files){
    names <- basename(i)
    name <- strsplit(names, split = "_")[[1]][1]
    colnames <- c(colnames, name)
}

input <- lapply(X = files, FUN = read.delim, header = FALSE, sep = "\t")
input.merge <- input %>% reduce(left_join, by = "V1")
colnames(input.merge) <- colnames

write.table(
    input.merge,
    file = outname,
    sep = "\t",
    quote = FALSE,
    col.names = NA,
    row.names = TRUE
)
message("combine counts complete!")
