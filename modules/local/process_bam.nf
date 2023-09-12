// Filters reads from BAM file that contain cell barcode (only contains R2 from CR)
// converts to fastq file
// https://kb.10xgenomics.com/hc/en-us/articles/360022448251-How-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
process PROCESS_BAM {
    tag "${sample_id}"
    label "process_high"
    tag "$sample_id"
    publishDir "${params.outdir}/process_bam/", mode: 'symlink'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}_R2.fastq.gz"), emit: reads

    script:
        """
        # Save the header lines
        samtools view -H $bam > SAM_header.sam

        # Filter alignments. Use LC_ALL=C to set C locale instead of UTF-8
        # can only filter for one tag in samtools
        # Combine header and body
        # convert BAM to fastq. CR output only contains R2
        # pipe everything to save time on IO
        cat SAM_header.sam  <(samtools view -d CB $bam | LC_ALL=C grep 'UB:Z:') |\
        samtools fastq -@ ${task.cpus} -0 ${sample_id}_R2.fastq.gz
        """
}
