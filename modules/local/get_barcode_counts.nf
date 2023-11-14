process GET_BARCODE_COUNTS {
    tag "$sample_id"
    label "process_low_bulk"

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        path "${sample_id}_rawcounts.txt"
    
    shell:
        """
        samtools idxstats -@ ${task.cpus} ${bam} | cut -f1,3 | awk '\$2!=0' > ${sample_id}_rawcounts.txt
        """
}