process BOWTIE_ALIGN {
    tag { "bowtie on ${sample_id}" }
    label "process_low"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
    tuple val(refname), path(ref_files)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.mapped.sam"), emit: mapped_reads
    path "${sample_id}.unmapped.fastq", emit: unmapped_reads
    path "${sample_id}.bowtie.log", emit: log

    script:
    """
    bowtie \\
    -p ${params.threads} \\
    -v ${params.alnmismatches} \\
    --norc \\
    -t \\
    ${refname} \\
    ${reads} \\
    --no-unal \\
    -S > ${sample_id}.mapped.sam \\
    2> ${sample_id}.bowtie.log
    """
}