// use cutadapt to filter for length
process CUTADAPT_READS{
  tag "cutadapt on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/trimmed_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq"), emit: reads
    path "${sample_id}.cutadapt.log", emit: log
  
  script:
  if(params.mode == "single-cell") {
    """
    cutadapt -j $((${params.threads} / 3)) -g "${params.upconstant}...${params.downconstant}" --trimmed-only --max-n=0  -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed_1.fastq 2> ${sample_id}.cutadapt_1.log &
    cutadapt -j $((${params.threads} / 3)) -g ${params.upconstant} --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed_2.fastq 2> ${sample_id}.cutadapt_2.log &
    cutadapt -j $((${params.threads} / 3)) -a ${params.downconstant} --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed_3.fastq 2> ${sample_id}.cutadapt_3.log &
    wait
    cat ${sample_id}.cutadapt_1.log ${sample_id}.cutadapt_2.log ${sample_id}.cutadapt_3.log > ${sample_id}.cutadapt.log
    cat ${sample_id}.trimmed_*.fastq > ${sample_id}.trimmed.fastq
    """
  }
  else if( params.constants == "both" )
    """
    cutadapt -j ${params.threads} -g "${params.upconstant}...${params.downconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
  else if( params.constants == "up" )
    """
    cutadapt -j ${params.threads} -g "${params.upconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
  else if( params.constants == "down" )
    """
    cutadapt -j ${params.threads} -a "${params.downconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
}