process COMBINE_STARCODE {
    tag { "combine starcode output files" }
    label "process_low"
    publishDir "${params.outdir}/counts/", mode: 'copy'

    input:
      path counts

    output:
      path "all_counts_combined.txt"
    
    script:
    """
    starcode_merge.py $counts all_counts_combined.txt
    """
}