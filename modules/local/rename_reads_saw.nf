process RENAME_READS_SAW {
    tag "Rename reads ${sample_id}"
    label "process_high"
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed_reads/", mode: 'symlink'

    // modify read name to keep MID in header
    // @DP8400029380TLL1C001R00701707112|Cx:i:37103|Cy:i:80544 57F58031B1F9 7B77C
    // to 
    // @DP8400029380TLL1C001R00701707112|37103_80544|7B77C

    input:
        // fastq with unmapped reads from SAW pipeline
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_renamed.fastq.gz")

    script:
        """
        gunzip -c ${reads} | parallel -j ${task.cpus} --pipe 'sed "s/Cx:i://g;s/|Cy:i:/_/g;s/ \\w* /|/g"' | gzip > ${sample_id}_renamed.fastq.gz
        """
}
