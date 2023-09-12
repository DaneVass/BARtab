// see here for syntax re: alignment indexes
// https://biocorecrg.github.io/SIB_course_nextflow_Nov_2021/docs/fourth.html
process BUILD_BOWTIE_INDEX {
    tag "$ref"
    label "process_medium"

    input:
        path ref

    output:
        tuple val("${ref}"), path("${ref}*.ebwt")

    script:
        """
        bowtie-build --threads ${task.cpus} ${ref} ${ref}
        """
}