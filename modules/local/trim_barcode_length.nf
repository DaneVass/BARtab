process TRIM_BARCODE_LENGTH{
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(sample_id), path(reads)

    output:
        // no need to publish these files, therefore no publishDir specified in modules.config
        tuple val(sample_id), path("${sample_id}.trimmed.${params.min_readlength}.fastq.gz"), emit: reads
    
    script:
        // Must write decompressed file, otherwise EOFError if piping gunzip output directly into fastx_trimmer
        """
            gunzip -c ${reads} > ${sample_id}.tmp.fastq

            if [ ${params.constant} == "up" ]; then
                # if upstream constant is trimmed, trim sequences on the 3' (right) end
                fastx_trimmer -z -l ${params.min_readlength} \\
                -i ${sample_id}.tmp.fastq \\
                -o ${sample_id}.trimmed.${params.min_readlength}.fastq.gz
            elif [ ${params.constant} == "down" ]; then
            # 2> ${sample_id}.trim.log \\
                # if downstream constant is trimmed, need to reverse reads to effectively trim the 5' end
                # write stderr to log file and pass stdout to next command
                fastx_reverse_complement -i ${sample_id}.tmp.fastq |\\
                fastx_trimmer -l ${params.min_readlength} |\\
                fastx_reverse_complement -z -o ${sample_id}.trimmed.${params.min_readlength}.fastq.gz
            fi         

            rm ${sample_id}.tmp.fastq
        """
}