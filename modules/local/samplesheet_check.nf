process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv

    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in bartab/bin/
    // params.input_type default is fastq
    // pass parameters if they have been set
    pipeline = params.pipeline == null ? "" : "-p ${params.pipeline}"
    // params.reference is whether or not to align to a reference. The reference for each sample is given in the sample sheet. 
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv \\
        ${params.mode} \\
        ${params.reference} \\
        -t ${params.input_type} \\
        $pipeline
    """
}
