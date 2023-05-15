#!/usr/bin/env Rscript

# INPUT
# - takes 2 command line arguments, 1. path to counts and 2. output path for tabulated barcode counts
# OUTPUT
# - a table with UMI counts per barcode per cell
# - 2 plots in pdf: barcodes_per_cell.pdf and UMIs_per_bc.pdf

library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)

############
##  Anno  ##
############

args <- commandArgs(TRUE)

dat <- args[1]
outfile <- args[2]

# read in barcode-cell pairs
dat <- read.delim(dat, header = TRUE)
dat <- dat[, c(2, 1, 3)]

# aggregate over barcodes and UMI counts using data.table
bc.counts <- as.data.frame(dcast(
  setDT(dat),
  cell ~ .,
  value.var = c("gene", "count"),
  fun.aggregate = function(x) {
    paste(x, collapse = ",")
  }
))

# Add "-1" to cell ID for metadata annotation purposes
bc.counts$cell <- paste0(bc.counts$cell, "-1")
rownames(bc.counts) <- bc.counts$cell
bc.counts$cell <- NULL
names(bc.counts) <- c("barcode", "bc.umi.count")

# write final output
write.table(bc.counts,
            outfile,
            quote = FALSE,
            row.names = TRUE)

###########
## Plots ##
###########

# Detected lineage barcodes per cell
lineagePerCell.dist.df <- dat %>%
  dplyr::select(cell, gene) %>%
  dplyr::group_by(cell) %>%
  dplyr::tally(., name = "number_of_lineage_barcodes") %>%
  dplyr::arrange(dplyr::desc(number_of_lineage_barcodes))

p <- ggplot(lineagePerCell.dist.df, aes(x = number_of_lineage_barcodes)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_histogram(
    aes(y = ..count.. / sum(..count..)),
    fill = "blue",
    binwidth = 0.5,
    alpha = 0.5
  ) +
  xlab("Detected barcodes per cell") +
  ylab("Fraction of cells") +
  ggtitle("Detected barcodes per cell") +
  theme_classic()

ggsave("barcodes_per_cell.pdf", p)

# Number of UMIs supporting the most frequent barcode
max.umi.per.cell <- setDT(dat)[, .(max = max(count)), by = list(cell)]

p <- ggplot(max.umi.per.cell, aes(x = max)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_histogram(aes(y = ..count.. / sum(..count..)), binwidth = 0.6) +
  xlab("UMI supporting the most frequent barcode") +
  ylab("Fraction of cells") +
  ggtitle("UMI supporting the most frequent barcode") +
  theme_classic()

ggsave("UMIs_per_bc.pdf", p)
