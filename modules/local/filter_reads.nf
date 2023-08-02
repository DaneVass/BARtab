process FILTER_READS{
  tag "fastx-toolkit fastq_filter_quality on $sample_id"
  label "process_medium"
  publishDir "${params.outdir}/filtered_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: reads
    path "${sample_id}.filter.log", emit: log
  
  script:
  // Must write decompressed file, otherwise EOFError if piping gunzip output directly into fastq_quality_filter
  """
  gunzip -c ${reads} > ${sample_id}.tmp.fastq
  fastq_quality_filter -i ${sample_id}.tmp.fastq \\
    -z -v -p ${params.pctqual} \\
    -q ${params.minqual} \\
    > ${sample_id}.filtered.fastq.gz \\
    2> ${sample_id}.filter.log

  rm ${sample_id}.tmp.fastq
  """
}