process GET_BARCODE_COUNTS {
    tag "$sample_id"
    label "process_low_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        path "${sample_id}_rawcounts.txt"
    
    shell:
        """
        samtools sort -@ ${task.cpus} ${reads} -o ${sample_id}.mapped_sorted.bam
        samtools index -@ ${task.cpus} ${sample_id}.mapped_sorted.bam ${sample_id}.mapped_sorted.bam.bai
        samtools idxstats -@ ${task.cpus} ${sample_id}.mapped_sorted.bam | cut -f1,3 | awk '\$2!=0' > ${sample_id}_rawcounts.txt
        """
}
