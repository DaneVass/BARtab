process FASTQC {
  tag "FastQC on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/qc/$sample_id", mode: 'copy'
  
  input:
    tuple val(sample_id), path(reads)

  output:
    tuple path("*.html"), path("*.zip")
  
  script:

  """
  fastqc --threads ${task.cpus} ${reads}
  """
}