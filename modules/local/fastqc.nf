process FASTQC {
    tag "$sample_id"
    label "$params.mode" == "single-cell" ? "process_low_sc" : "process_low_bulk"
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple path("${sample_id}*.html"), path("${sample_id}*.zip")
    
    script:
        """
        fastqc --threads ${task.cpus} ${reads}
        """
}