process MERGE_READS {
  tag "FLASh on $sample_id"
  label "process_high"
  publishDir "${params.outdir}/merged_reads/$sample_id", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: merged_reads
    path "${sample_id}.notCombined_1.fastq.gz"
    path "${sample_id}.notCombined_2.fastq.gz"
    path "${sample_id}.hist"
    path "${sample_id}.histogram"
    path "${sample_id}.flash.log", emit: log

  script: 
  """
  flash -z --min-overlap=${params.mergeoverlap} -t ${params.threads} --output-prefix=${sample_id} ${reads} 2> "${sample_id}.flash.log"
  """
}
