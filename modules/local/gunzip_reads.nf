process GUNZIP_READS{
  tag "gunzip -f on $sample_id"
  label "process_low"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*.fastq")
  
  script:
  """
  gunzip -f ${reads}
  """
}