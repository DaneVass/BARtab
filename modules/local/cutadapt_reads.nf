// use cutadapt to filter for length
process CUTADAPT_READS{
    tag "$sample_id"
    label "$params.mode" == "single-cell" ? "process_high_sc" : "process_high_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        // possibility that no reads contain adapter, will be handled by not zipping empty output file. This way, pipeline does not crash. 
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: reads
        path "${sample_id}.cutadapt_*.log", emit: log
    
    script:

        // filter for max length if barcode has fixed length
        def max_len = params.barcode_length ? "-M ${params.barcode_length}" : ""
        // filter for whole barcode length if both adapters present if barcode has fixed length
        // bowtie performs ungapped alignment
        def min_len_both = params.barcode_length ? "-m ${params.barcode_length}" : "-m ${params.min_readlength}"

        if(params.constants == "all") {
            """
            cutadapt -j \$((${task.cpus} / 3)) -g "${params.upconstant};o=${params.up_coverage};required...${params.downconstant};o=${params.down_coverage};required" --trimmed-only --max-n=0  ${min_len_both} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_both.fastq 2> ${sample_id}.cutadapt_both.log &
            cutadapt -j \$((${task.cpus} / 3)) -g "${params.upconstant};o=${params.up_coverage}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_up.fastq 2> ${sample_id}.cutadapt_up.log &
            cutadapt -j \$((${task.cpus} / 3)) -a "${params.downconstant};o=${params.down_coverage}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed_down.fastq 2> ${sample_id}.cutadapt_down.log &
            wait
            cat ${sample_id}.trimmed_*.fastq > ${sample_id}.trimmed.fastq

            # check if cutadapt trimmed any reads, else throw error to exclude file
            # if not, bowtie would crash with a hard to interpret error message
            if [ -s ${sample_id}.trimmed.fastq ]; then
                pigz -p ${task.cpus} ${sample_id}.trimmed.fastq
            else
                printf '%s\n' "No reads passed cutadapt filtering and trimming. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
                exit 1
            fi
            """
        }
        else if( params.constants == "both" )
            """
            cutadapt -j ${task.cpus} -g "${params.upconstant};o=${params.up_coverage};required...${params.downconstant};o=${params.down_coverage};required" --trimmed-only --max-n=0 ${min_len_both} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_both.log
            if [ -s ${sample_id}.trimmed.fastq ]; then
                pigz -p ${task.cpus} ${sample_id}.trimmed.fastq
            else
                printf '%s\n' "No reads passed cutadapt filtering and trimming. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
                exit 1
            fi
            """
        else if( params.constants == "up" )
            """
            cutadapt -j ${task.cpus} -g "${params.upconstant};o=${params.up_coverage}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_up.log
            if [ -s ${sample_id}.trimmed.fastq ]; then
                pigz -p ${task.cpus} ${sample_id}.trimmed.fastq
            else
                printf '%s\n' "No reads passed cutadapt filtering and trimming. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
                exit 1
            fi
            """
        else if( params.constants == "down" )
            """
            cutadapt -j ${task.cpus} -a "${params.downconstant};o=${params.down_coverage}" --trimmed-only --max-n=0 -m ${params.min_readlength} -e ${params.constantmismatches} ${max_len} ${reads} > ${sample_id}.trimmed.fastq 2> ${sample_id}.cutadapt_down.log
            if [ -s ${sample_id}.trimmed.fastq ]; then
                pigz -p ${task.cpus} ${sample_id}.trimmed.fastq
            else
                printf '%s\n' "No reads passed cutadapt filtering and trimming. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
                exit 1
            fi
            """
}