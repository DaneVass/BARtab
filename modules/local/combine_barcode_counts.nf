process COMBINE_BARCODE_COUNTS {
    label "process_medium_bulk"

    input:
      path counts

    output:
      path "all_counts_combined.tsv"

    script:
      """
      combine_counts.py $counts all_counts_combined.tsv
      """
}