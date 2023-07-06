process RENAME_READS_SAW {
    tag { "Rename reads ${sample_id}" }
    label "process_high"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    // modify read name so umi_tools can handle it
    // from this DP8400022151TRL1C001R00400077156|||CB:Z:48489_71597|||UR:Z:CGCTTGGCCT|||UY:Z:FFEECFEEEG
    // to this DP8400022151TRL1C001R00400077156|||48489_71597|||CGCTTGGCCT

    input:
    // sam containing mapped reads
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.mapped_renamed.bam")

    script:
    """
    parallel -j ${task.cpus} -a ${sam} --pipepart 'sed "s/CB:Z://g;s/UR:Z://g;s/|||UY:Z:\\S*//g"' |\
      samtools view -@ ${task.cpus} -b -o ${sample_id}.mapped_renamed.bam
    """
}
