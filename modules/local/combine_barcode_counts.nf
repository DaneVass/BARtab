process COMBINE_BARCODE_COUNTS {
  label "process_low"
  publishDir "${params.outdir}/counts/", mode: 'copy'

  input:
    path counts

  output:
    path "all_counts_combined.txt"
  
  script: 
  """
  combine_counts.R $counts all_counts_combined.txt
  """
}