process GUNZIP_READS_PE {
  tag "gunzip -f on $sample_id"
  label "process_low"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq")
  
  script:
  """
  gunzip -f ${reads}
  """
}