// add cell barcode and UMI to read name
// performed at this point only on aligned sequences because the process is very slow on large files
process RENAME_READS {
    tag "Rename reads ${sample_id}"
    label "process_high"
    tag "$sample_id"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
        // sam containing mapped reads, bam from cellranger
        tuple val(sample_id), path(sam), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.mapped_renamed.sam")

    script:
        """
        # Save the header lines
        samtools view -H $sam > SAM_header.sam

        # get lines from input bam (containing cb and umi info) that are mapped to ref
        # create header that contains cb and umi, needed for umi_tools count
        # first column is original header to merge on later
        # piping awk output into sed is faster and uses less CPU and IO and only marginally more memory
        awk -F"\t" 'FNR==NR {lines[\$1]; next} \$1 in lines' <(samtools view $sam) <(samtools view $bam) |\
        parallel -j ${task.cpus} --pipe 'sed "s/^\\([^\t]*\\).*CB:Z:\\([^\t]*\\)-1.*UB:Z:\\([^\t]*\\).*/\\1\t\\1_\\2_\\3/g"' |\
        sort -S 2G -k1,1 > read_names.txt

        # join based on first column, then print all but first column
        join -j1 -t "\$(echo -e "\t")" read_names.txt <(samtools view $sam | sort -S 2G -k1,1) |\\
        awk -v OFS="\t" '{\$1=""; print substr(\$0,2)}' > SAM_body_renamed.sam

        # Combine header and body, compress
        cat SAM_header.sam SAM_body_renamed.sam > ${sample_id}.mapped_renamed.sam

        rm SAM_header.sam read_names.txt SAM_body_renamed.sam
        """
}
