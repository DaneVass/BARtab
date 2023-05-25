process UMITOOLS_COUNT {
  tag "umi_tools count on $sample_id"
  label "process_medium"
  publishDir "${params.outdir}/counts", mode: 'symlink'

  input:
    tuple val(sample_id), path(bam)
    // index file needs to be linked to work directory
    path bai

  output:
    tuple val(sample_id), path("${sample_id}.counts.tsv")
    
  script:
    """
    umi_tools count \\
    --per-contig --per-cell \\
    --edit-distance-threshold=${params.umi_dist} \\
    --random-seed=10101 \\
    -I ${bam} \\
    -S ${sample_id}.counts.tsv
    """  
}