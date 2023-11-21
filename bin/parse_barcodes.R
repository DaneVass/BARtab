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
umi_fraction_filter <- args[4]
umi_count_filter <- args[5]

umi_fraction_filter <- as.numeric(umi_fraction_filter)
umi_count_filter <- as.numeric(umi_count_filter)

# read in barcode-cell pairs
counts_data <- read.delim(counts_data, header = TRUE)
counts_data <- counts_data[, c(2, 1, 3)]

counts_data_filtered <- counts_data %>%
  filter(count >= umi_count_filter) %>%
  group_by(cell) %>%
  mutate(max_count = max(count)) %>%
  arrange(cell) %>%
  filter(count / umi_fraction_filter >= max_count) %>%
  dplyr::select(-max_count) %>%
  ungroup()

# aggregate over barcodes and UMI counts using data.table
bc.counts <- as.data.frame(dcast(
  setDT(counts_data_filtered),
  cell ~ .,
  value.var = c("gene", "count"),
  fun.aggregate = function(x) {
    paste(x, collapse = ";")
  }
))

rownames(bc.counts) <- bc.counts$cell
bc.counts$cell <- NULL
names(bc.counts) <- c("barcode", "bc.umi.count")

# write final output
write.table(
  bc.counts,
  paste0(sample_id, "_cell_barcode_annotation.tsv"),
  quote = FALSE,
  row.names = TRUE,
  sep = "\t"
)

###########
## Plots ##
###########

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

# Detected lineage barcodes per cell
lineagePerCell.dist.df <- counts_data %>%
  dplyr::select(cell, gene) %>%
  dplyr::group_by(cell) %>%
  dplyr::tally(., name = "number_of_lineage_barcodes") %>%
  dplyr::arrange(dplyr::desc(number_of_lineage_barcodes))

p <- lineagePerCell.dist.df %>%
  dplyr::count(number_of_lineage_barcodes) %>%
  mutate(frac = n / nrow(lineagePerCell.dist.df)) %>%
  ggplot(aes(x = number_of_lineage_barcodes, y = frac)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_bar(stat = "identity",
           fill = "blue",
           width = .5) +
  xlab("Detected barcodes per cell") +
  ylab("Fraction of cells") +
  ggtitle("Detected barcodes per cell") +
  theme_classic() +
  scale_x_continuous(breaks = integer_breaks())

ggsave(paste0(sample_id, "_barcodes_per_cell.pdf"), p)

# Detected lineage barcodes per cell
lineagePerCell.dist.df.filtered <- counts_data_filtered %>%
  dplyr::select(cell, gene) %>%
  dplyr::group_by(cell) %>%
  dplyr::tally(., name = "number_of_lineage_barcodes") %>%
  dplyr::arrange(dplyr::desc(number_of_lineage_barcodes))

p <- lineagePerCell.dist.df.filtered %>%
  dplyr::count(number_of_lineage_barcodes) %>%
  mutate(frac = n / nrow(lineagePerCell.dist.df)) %>%
  ggplot(aes(x = number_of_lineage_barcodes, y = frac)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_bar(stat = "identity",
           fill = "blue",
           width = .5) +
  xlab("Detected barcodes per cell") +
  ylab("Fraction of cells") +
  ggtitle("Detected barcodes per cell, filtered") +
  theme_classic() +
  scale_x_continuous(breaks = integer_breaks())

ggsave(paste0(sample_id, "_barcodes_per_cell_filtered.pdf"), p)

# Number of UMIs supporting the most frequent barcode
max.umi.per.cell <-
  setDT(counts_data)[, .(max = max(count)), by = list(cell)]

p <- max.umi.per.cell %>%
  dplyr::count(max) %>%
  mutate(frac = n / nrow(max.umi.per.cell)) %>%
  ggplot(aes(x = max, y = frac)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_bar(stat = "identity") +
  xlab("UMI supporting the most frequent barcode") +
  ylab("Fraction of cells") +
  ggtitle("UMI supporting the most frequent barcode") +
  coord_cartesian(xlim = c(0, NA)) +
  theme_classic() +
  scale_x_continuous(breaks = integer_breaks())

ggsave(paste0(sample_id, "_UMIs_per_bc.pdf"), p)

max.umi.per.cell.filtered <-
  setDT(counts_data_filtered)[, .(max = max(count)), by = list(cell)]

p <- max.umi.per.cell.filtered %>%
  dplyr::count(max) %>%
  mutate(frac = n / nrow(max.umi.per.cell)) %>%
  ggplot(aes(x = max, y = frac)) +
  theme(
    axis.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold")
  ) +
  geom_bar(stat = "identity") +
  xlab("UMI supporting the most frequent barcode") +
  ylab("Fraction of cells") +
  ggtitle("UMI supporting the most frequent barcode, filtered") +
  coord_cartesian(xlim = c(0, NA)) +
  theme_classic() +
  scale_x_continuous(breaks = integer_breaks())

ggsave(paste0(sample_id, "_UMIs_per_bc_filtered.pdf"), p)

if (file.size(sam_file) != 0L) {
  # read SAM file of aligned sequences
  sam <-
    read.delim(sam_file,
              sep = "\t",
              comment.char = "@",
              header = F)
  # select columns of read name, barcode and sequence
  sam <- unique(sam[, c(1, 3, 10)])
  sam <- sam[, c(2, 3)]
  colnames(sam) <- c("barcode", "sequence")
  # calculate average sequence length per barcode
  avg_seq_len <- sam %>%
    mutate(sequence_length = nchar(as.character(sequence))) %>%
    select(-sequence) %>%
    group_by(barcode) %>%
    summarize(avg_sequence_length = mean(sequence_length)) %>%
    arrange(desc(avg_sequence_length))

  # save data to be able to filter barcodes on their avg seq length later on

  write.table(
    avg_seq_len,
    paste0(sample_id, "_avg_sequence_length.tsv"),
    quote = FALSE,
    row.names = TRUE,
    sep = "\t"
  )

  # count in how many cells each barcode is detected, allows for detection of multiple barcodes in a cell
  counts_agg <-
    data.frame(table(counts_data$gene))
  colnames(counts_agg) <- c("barcode", "count")

  # merge count and sequence length data and plot
  p <- merge(counts_agg, avg_seq_len, by = "barcode") %>%
    ggplot(aes(x = count, y = avg_sequence_length)) +
    geom_point(alpha = 0.2) +
    ggtitle("Average length of sequences mapped to barcode reference") +
    ylab("Average mapped sequence length") +
    xlab("cells")

  # save plot
  ggsave(paste0(sample_id, "_avg_sequence_length.pdf"), p)

  counts_agg <-
    data.frame(table(counts_data_filtered$gene))
  colnames(counts_agg) <- c("barcode", "count")

  # merge count and sequence length data and plot
  p <- merge(counts_agg, avg_seq_len, by = "barcode") %>%
    ggplot(aes(x = count, y = avg_sequence_length)) +
    geom_point(alpha = 0.2) +
    ggtitle("Average length of sequences mapped to barcode reference, filtered") +
    ylab("Average mapped sequence length") +
    xlab("cells")

  # save plot
  ggsave(paste0(sample_id, "_avg_sequence_length_filtered.pdf"), p)
}