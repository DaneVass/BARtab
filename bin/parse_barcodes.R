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

counts_data <- args[1]
sample_id <- args[2]
sam_file <- args[3]

# read in barcode-cell pairs
counts_data <- read.delim(counts_data, header = TRUE)
counts_data <- counts_data[, c(2, 1, 3)]

# aggregate over barcodes and UMI counts using data.table
bc.counts <- as.data.frame(dcast(
  setDT(counts_data),
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
            paste0(sample_id, "_cell-barcode-anno.tsv"),
            quote = FALSE,
            row.names = TRUE)

###########
## Plots ##
###########

# Detected lineage barcodes per cell
lineagePerCell.dist.df <- counts_data %>%
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

ggsave(paste0(sample_id, "_barcodes_per_cell.pdf"), p)

# Number of UMIs supporting the most frequent barcode
max.umi.per.cell <- setDT(counts_data)[, .(max = max(count)), by = list(cell)]

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

ggsave(paste0(sample_id, "_UMIs_per_bc.pdf"), p)


# read SAM file of aligned sequences
sam <- read.delim(sam_file, sep ="\t", comment.char = "@", header = F)
# select columns of read name, barcode and sequence
sam <- unique(sam[,c(1,3,10)])
sam <- sam[,c(2,3)]
colnames(sam) <- c("barcode", "sequence")
# calculate average sequence length per barcode
avg_seq_len <- sam %>%
  mutate(sequence_length = nchar(as.character(sequence))) %>%
  select(-sequence) %>%
  group_by(barcode) %>%
  summarize(avg_sequence_length = mean(sequence_length)) %>%
  arrange(desc(avg_sequence_length))

# save data to be able to filter barcodes on their avg seq length later on

write.table(avg_seq_len,
            paste0(sample_id, "_avg_sequence_length.tsv"),
            quote = FALSE,
            row.names = TRUE)

# count in how many cells each barcode is detected, allows for detection of multiple barcodes in a cell
counts_agg <-
  data.frame(table(counts_data$gene))
colnames(counts_agg) <- c("barcode", "count")

# merge count and sequence length data and plot
p <- merge(counts_agg, avg_seq_len, by="barcode") %>%
  ggplot(aes(x=count, y=avg_sequence_length)) +
  geom_point(alpha=0.2) +
  ggtitle("Average length of sequences mapped to barcode reference") +
  ylab("Average mapped sequence length") +
  xlab("UMI count")

# save plot
ggsave(paste0(sample_id, "_avg_sequence_length.pdf"), p)
