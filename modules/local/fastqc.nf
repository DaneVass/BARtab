process FASTQC {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/qc/$sample_id", mode: 'copy'
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple path("${sample_id}*.html"), path("${sample_id}*.zip")
    
    script:
        """
        fastqc --threads ${task.cpus} ${reads}
        """
}