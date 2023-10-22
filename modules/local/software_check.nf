// check and report on software versions used in the pipeline run
process SOFTWARE_CHECK {
    label 'process_single'

    output:
        path("software_check.txt")

    script:
        """
        check_versions.sh software_check.txt
        """
}