process GET_BARCODE_COUNTS {
  tag "samtools idxstats on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/counts/", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    path "${sample_id}_rawcounts.txt"
  
  shell:
  """
  samtools idxstats -@ ${task.cpus} ${reads} | cut -f1,3 > ${sample_id}_rawcounts.txt
  """
}