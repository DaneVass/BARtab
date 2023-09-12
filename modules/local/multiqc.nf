process MULTIQC {
    label "process_low"

    publishDir "${params.outdir}", mode: 'copy', overwrite: 'true'

    input:
        path multiqcConfig
        path output
        path counts

    output:
        path "multiqc_report.html"
        path "multiqc_data"

    script:
        """
        multiqc $output -f
        """
}