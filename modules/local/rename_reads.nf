// add cell barcode and UMI to read name
// performed at this point only on aligned sequences because the process is very slow on large files
process RENAME_READS {
    tag { "Rename reads ${sample_id}" }
    label "process_high"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
    // sam containing mapped reads
    tuple val(sample_id), path(sam)
    // bam from cellranger
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.mapped_renamed.bam")

    script:
    """
    # Save the header lines
    samtools view -H $sam > SAM_header.sam

    # get lines from input bam (containing cb and umi info) that are mapped to ref
    # takes 5min on 8.5GB bam
    awk -F"\t" 'FNR==NR {lines[\$1]; next} \$1 in lines' <(samtools view $sam) <(samtools view $bam) > ${sample_id}_rename_helper.sam

    # create header that contains cb and umi, needed for umi_tools count
    # first column is original header to merge on later
    # parallelize sed command, otherwise takes 15min on 10,450,943 mapped reads, with parallel and 16 cores only 2min
    parallel -j ${task.cpus} -a ${sample_id}_rename_helper.sam --pipepart 'sed "s/^\\([^\t]*\\).*CB:Z:\\([^\t]*\\)-1.*UB:Z:\\([^\t]*\\).*/\\1\t\\1_\\2_\\3/g"' | sort -S 2G -k1,1 > read_names.txt

    # join based on first column, then print all but first column
    join -j1 -t "\$(echo -e "\t")" read_names.txt <(samtools view $sam | sort -S 2G -k1,1) |\\
      awk -v OFS="\t" '{\$1=""; print substr(\$0,2)}' > SAM_body_renamed.sam

    # Combine header and body
    cat SAM_header.sam SAM_body_renamed.sam | samtools view -@ ${task.cpus} -b -o ${sample_id}.mapped_renamed.bam

    rm SAM_header.sam ${sample_id}_rename_helper.sam read_names.txt SAM_body_renamed.sam
    """
}
