process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    // conda "conda-forge::python=3.8.3"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/python:3.8.3' :
    //     'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv

    when:
    task.ext.when == null || task.ext.when

    // pass parameters if they have been set
    pipeline = params.pipeline == null ? "" : "-p ${params.pipeline}"
    input_type = params.input_type == null ? "" : "-t ${params.input_type}"

    script: // This script is bundled with the pipeline, in bartab/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv \\
        ${params.mode} \\
        ${params.reference} \\
        $input_type \\
        $pipeline
    """
}
