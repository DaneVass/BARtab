

// SINGLE CELL WORKFLOW BEGINS HERE
// 05A_umitools whitelist && extract
if( params.singlecell) {
    println "Begin single-cell workflow"
    println ""
    Channel
      .fromPath( "${params.indir}/*_R1.{fastq,fq}.gz" )
      .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
      .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }
      .into { scChannel }
}

process umitools_whitelist{
  tag "umi_tools whitelist on $sample_id"
  label "process_med"
  publishDir "${params.outdir}/umi_tools", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads) from scChannel

  output:
    set val(sample_id), file("${sample_id}_whitelist.txt") into scWhitelistChannel
    file("${sample_id}.cutadapt.log") into trimmedLogChannel

  script:
  if( params.singlecell )
    """
    umi_tools whitelist --stdin ${reads} \\
      --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \\
      --log2stderr > ${sample_id}_whitelist.txt
    """  
}

process umitools_extract{
  tag "umi_tools extract on $sample_id"
  label "process_med"
  publishDir "${params.outdir}/umi_tools", mode: 'symlink'

  input:
    tuple val(sample_id), path(whitelist) from scWhitelistChannel
    tuple val(sample_id), path(reads) from filteredReadsChannel

  output:
    set val(sample_id), file("${sample_id}_R2_extracted.fastq") into scExtractChannel
    file("${sample_id}.cutadapt.log") into filteredReadsChannel

  script:
  if( params.singlecell )
    """
    umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \\
      --stdin ${reads[0]} --stdout extract/${sample_id}_extracted.fastq \\
      --read2-in ${reads[1]} --read2-out=extract/${sample_id}_R2_extracted.fastq \\
      --whitelist=${whitelist}
    """  
}