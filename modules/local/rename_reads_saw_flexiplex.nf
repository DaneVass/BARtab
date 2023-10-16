process RENAME_READS_FLEXIPLEX {
    tag "Rename reads ${sample_id}"
    label "process_high"
    tag "$sample_id"

    // modify read name from FLAMES/flexiplex output
    // @CGTAAGTCAGTAGGAC_TCAAGTAGGTGG#147deabb-756c-48ea-a88c-571514f4abdb_-1of1
    // to 
    // @147deabb-756c-48ea-a88c-571514f4abdb_-1of1_CGTAAGTCAGTAGGAC_TCAAGTAGGTGG

    input:
        // fastq with unmapped reads from SAW pipeline
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_renamed.fastq.gz")

    script:
        """
        gunzip -c ${reads} | parallel -j ${task.cpus} --pipe 'sed "'s/@\\([ACTG_]*\\)#\\(.*\\)/@\\2_\\1/g'"' | gzip > ${sample_id}_renamed.fastq.gz
        """
}
