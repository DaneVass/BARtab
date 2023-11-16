process TRIM_BARCODE_LENGTH{
    tag "$sample_id"
    label "process_medium_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        // no need to publish these files, therefore no publishDir specified in modules.config
        tuple val(sample_id), path("${sample_id}.trimmed.${params.min_readlength}.fastq.gz"), emit: reads
    
    script:
        // Must write decompressed file, otherwise EOFError if piping gunzip output directly into fastx_trimmer
        """
        if [ ${params.constants} == "up" ]; then
            # if upstream constant is trimmed, trim sequences on the 3' (right) end
            fastp --in1 ${reads} \\
                --thread ${task.cpus} \\
                --out1 ${sample_id}.trimmed.${params.min_readlength}.fastq.gz \\
                --disable_trim_poly_g --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --dont_eval_duplication \\
                --max_len1 ${params.min_readlength}
        elif [ ${params.constants} == "down" ]; then
            # if downstream constant is trimmed, need to reverse reads to effectively trim the 5' end
            # then reverse reads again
            # I haven't yet encountered a case where only the downstream adapter is trimmed and barcodes must be trimmed for clustering
            # there might be a more elegant way to do this
            zcat ${reads} |\\
                awk '!(FNR%2) {t="";l=length;for(i=1;i<=l;i++) t=substr(\$0,i,1) t;\$0=t}1' |\\
                fastp --stdin --stdout \\
                    --thread ${task.cpus} \\
                    --disable_trim_poly_g --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --dont_eval_duplication \\
                    --max_len1 ${params.min_readlength} |\\
                awk '!(FNR%2) {t="";l=length;for(i=1;i<=l;i++) t=substr(\$0,i,1) t;\$0=t}1' |\\
                gzip > ${sample_id}.trimmed.${params.min_readlength}.fastq.gz
        fi
        """
}