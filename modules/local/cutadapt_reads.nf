// use cutadapt to filter for length
process CUTADAPT_READS{
    tag "$sample_id"
    label "process_high"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: reads
        path "${sample_id}.cutadapt_*.log", emit: log
    
    script:

        // filter for max length if barcode has fixed length
        def max_len = params.barcode_length ? "-M ${params.barcode_length}" : ""
        // filter for whole barcode length if both adapters present if barcode has fixed length
        // bowtie performs ungapped alignment
        def min_len_both = params.barcode_length ? "-m ${params.barcode_length}" : "-m ${params.min_readlength}"

        if(params.mode == "single-cell" || params.constants == "all") {
            """
            cutadapt -j \$((${task.cpus} / 3)) -g "${params.upconstant};required...${params.downconstant};required" --trimmed-only --max-n=0  ${min_len_both} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_1.fastq 2> ${sample_id}.cutadapt_both.log &
            cutadapt -j \$((${task.cpus} / 3)) -g "${params.upconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_2.fastq 2> ${sample_id}.cutadapt_up.log &
            cutadapt -j \$((${task.cpus} / 3)) -a "${params.downconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_3.fastq 2> ${sample_id}.cutadapt_down.log &
            wait
            cat ${sample_id}.trimmed_*.fastq > ${sample_id}.trimmed.fastq
            gzip ${sample_id}.trimmed.fastq
            """
        }
        else if( params.constants == "both" )
            """
            cutadapt -j ${task.cpus} -g "${params.upconstant};required...${params.downconstant};required" --trimmed-only --max-n=0 ${min_len_both} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_both.log
            gzip ${sample_id}.trimmed.fastq
            """
        else if( params.constants == "up" )
            """
            cutadapt -j ${task.cpus} -g "${params.upconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_up.log
            gzip ${sample_id}.trimmed.fastq
            """
        else if( params.constants == "down" )
            """
            cutadapt -j ${task.cpus} -a "${params.downconstant}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_down.log
            gzip ${sample_id}.trimmed.fastq
            """
}