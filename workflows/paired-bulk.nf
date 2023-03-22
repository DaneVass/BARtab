// BARtab paired end bulk workflow

//--------------------------------------------------------------------------------------
// Main Pipeline 
//--------------------------------------------------------------------------------------

//process 1 file unless --n used at run time, e.g. -n 32
params.n = 1 

// 00 Process inputs - assume paired-end reads
Channel
  .fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
  .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.reads}" }
  .set { readsChannel }

readsChannel.view { "file: $it" }

// 01_FastQC
process fastqc{
  tag "FastQC on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/qc/$sample_id", mode: 'copy'
  
  input:
    tuple val(sample_id), path(reads)

  output:
    tuple file('*.html'), file('*.zip')
  
  script:

  """
  fastqc --threads ${params.threads} ${reads}
  """
}

//02_Read_merging
process merge_reads{
  tag "FLASh on $sample_id"
  label "process_high"
  publishDir "${params.outdir}/merged_reads/$sample_id", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), file("${sample_id}.extendedFrags.fastq.gz")
    file "${sample_id}.notCombined_1.fastq.gz"
    file "${sample_id}.notCombined_2.fastq.gz"
    file "${sample_id}.hist"
    file "${sample_id}.histogram"
    file "${sample_id}.flash.log"

  script: 
  """
  flash -z --min-overlap=10 -t ${params.threads} --output-prefix=${sample_id} ${reads} 2> "${sample_id}.flash.log"
  """
}

// 03_gunzip_files
process gunzip_reads{
  tag "gunzip -f on $sample_id"
  label "process_low"
  
  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), file("${sample_id}.extendedFrags.fastq")
  
  script:
  """
  gunzip -f ${reads}
  """
}

// 04_Filter_reads_on_quality
process filter_reads{
  tag "fastx-toolkit fastq_filter_quality on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/filtered_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), file("${sample_id}.filtered.fastq.gz")
    file("${sample_id}.filter.log")

  script:
  """
  fastq_quality_filter -z -v -p ${params.pctqual} -q ${params.minqual} -i ${reads} > ${sample_id}.filtered.fastq.gz 2> ${sample_id}.filter.log
  """
}

// 05_cutadapt // use cutadapt to filter for length
process cutadapt_reads{
  tag "cutadapt on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/trimmed_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), file("${sample_id}.trimmed.fastq")
    file("${sample_id}.cutadapt.log")
  
  script:
  if( params.merge )
    """
    cutadapt -g "${params.upconstant}...${params.downconstant}" --trimmed-only --max-n=0 -m 15 ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
  else if( params.constants == "both" )
    """
    cutadapt -g "${params.upconstant}...${params.downconstant}" --trimmed-only --max-n=0 -m 15 ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
  else if( params.merge && params.constants == "up" )
    """
    cutadapt -g "${params.upconstant}" --trimmed-only --max-n=0 -m 15 ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
  else
    """
    cutadapt -g "${params.upconstant}" --trimmed-only --max-n=0 -m 15 ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt.log
    """
}

// see here for syntax re: alignment indexes
// https://biocorecrg.github.io/SIB_course_nextflow_Nov_2021/docs/fourth.html

// 06_generate_bowtie_index
reference = file(params.ref)
process build_index {
    tag { "bowtie_build on ${ref}" }

    input:
    path ref

    output:
    tuple val("${ref}"), path ("${ref}*.ebwt")

    script:
    """
    bowtie-build ${ref} ${ref}
    """
}

/// 07_bowtie_align
process bowtie_align {
    tag { "bowtie on ${sample_id}" }
    label "process_low"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
    tuple val(refname), path (ref_files)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}.mapped.sam")
    file("${sample_id}.unmapped.fastq")
    file("${sample_id}.bowtie.log")

    script:
    """
    bowtie \\
    -p ${params.threads} \\
    -v ${params.alnmismatches} \\
    --norc \\
    -t \\
    ${refname} \\
    ${reads} \\
    --un ${sample_id}.unmapped.fastq \\
    -S > ${sample_id}.mapped.sam \\
    2> ${sample_id}.bowtie.log \\
    """
}

// 08_samtools
process samtools {
    tag { "samtools on ${sample_id}" }
    label "process_low"
    publishDir "${params.outdir}/mapped_reads/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}.mapped.bam")
    file("${sample_id}.mapped.bam.bai")
    
    script:
    """
    samtools view -Sb ${reads} | samtools sort - > ${sample_id}.mapped.bam
    samtools index ${sample_id}.mapped.bam ${sample_id}.mapped.bam.bai
    """
}

// 09_get_barcode_counts
process get_barcode_counts{
  tag "samtools idxstats on $sample_id"
  label "process_low"
  publishDir "${params.outdir}/counts/", mode: 'copy'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), file("${sample_id}_rawcounts.txt")
  
  shell:
  """
  samtools idxstats ${reads} | cut -f1,3 > ${sample_id}_rawcounts.txt
  """
}

// 10_combine_barcode_counts
process combine_barcode_counts{
  label "process_low"
  publishDir "${params.outdir}/counts/", mode: 'copy'

  input:
  file(counts)

  output:
    file "all_counts_combined.txt"
  
  script: 
  """
  Rscript $projectDir/scripts/combine_counts.R $counts all_counts_combined.txt
  """
}

// 11_multiqc_report
params.multiqc_config = "$projectDir/config/multiqc_config.yaml"
Channel
  .fromPath(params.multiqc_config, checkIfExists: true)
  .set { multiqcConfig }

process multiqc {
  label "process_low"

  publishDir "${params.outdir}", mode: 'copy', overwrite: 'true'

  input:
    path multiqcConfig
    file ("qc/*")
    file (fastx)
    file (cutadapt)
    file (bowtie)

  output:
    file "multiqc_report.html"
    file "multiqc_data"

  script:
  """
  multiqc . -f
  """

}

// 12_software_check
// check and report on software versions used in the pipeline run
process software_check {
  label 'software_check'

  publishDir params.outdir, mode: 'copy', overwrite: 'true'

  output:
    file("software_check.txt")

  script:
  """
  bash $projectDir/scripts/check_versions.sh software_check.txt
  """
}

//--------------------------------------------------------------------------------------
// Post processing
//--------------------------------------------------------------------------------------

// Mail notification

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info '\n'
}
else {
    log.info "\n"
    log.info "Sending runtime report to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()
    
    sendMail(to: params.email, subject: "BARtab execution report", body: msg,  attach: "${params.outdir}/multiqc_report.html")
    }
}

// Print completion messages
workflow.onComplete {
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'
    
    log.info " ---------------------- BARtab Pipeline has finished ----------------------"
    log.info ""
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Pipeline runtime: ${workflow.duration}\n"
}
