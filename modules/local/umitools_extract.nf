process UMITOOLS_EXTRACT {
  tag "umi_tools extract on $sample_id"
  label "process_medium"
  publishDir "${params.outdir}/extract", mode: 'symlink'

  input:
    tuple val(sample_id), path(whitelist), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_R2_extracted.fastq"), emit: reads
    path("${sample_id}_extract.log"), emit: log

  script:
    """
    umi_tools extract --bc-pattern=${params.cb_umi_pattern} \\
      --stdin ${reads[0]} \\
      --read2-in ${reads[1]} --read2-out=${sample_id}_R2_extracted.fastq \\
      --whitelist=${whitelist} -L ${sample_id}_extract.log
    """
}