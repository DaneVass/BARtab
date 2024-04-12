// add cell barcode and UMI to read name
// necessary for BAM from cellranger or STARsolo as input
// performed at this point only on aligned sequences because the process is very slow on large files
process RENAME_READS_SPLITPIPE {
    label "process_high_sc"
    tag "$sample_id"

    input:
        // sam containing mapped reads, bam from cellranger
        tuple val(sample_id), path(sam), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.mapped_renamed.sam")

    script:
        """
        # Save the header lines
        # samtools view -H $sam > SAM_header.sam

        cat <(samtools view -H $sam) \\
            <(samtools view $sam |\\
               sed 's/^\\([0-9]*_[0-9]*_[0-9]*\\)\\(.*\\)__\\([ACTG]*\\)__\\([^\t]*\\)/\\1\\2__\\3__\\4|\\1|\\3/g') \\
            > ${sample_id}.mapped_renamed.sam

        # Combine header and body, compress
        # cat SAM_header.sam SAM_body_renamed.sam > ${sample_id}.mapped_renamed.sam

        # rm SAM_header.sam read_names.txt SAM_body_renamed.sam
        """
}
