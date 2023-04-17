// see here for syntax re: alignment indexes
// https://biocorecrg.github.io/SIB_course_nextflow_Nov_2021/docs/fourth.html
// reference = file(params.ref)
process BUILD_BOWTIE_INDEX {
    tag { "bowtie_build on ${ref}" }

    input:
    path ref

    output:
    tuple val("${ref}"), path("${ref}*.ebwt")

    script:
    """
    bowtie-build ${ref} ${ref}
    """
}