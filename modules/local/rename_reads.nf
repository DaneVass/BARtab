// Filters reads from BAM file that contain cell barcode (only contains R2 from CR)
// converts to fastq file
// https://kb.10xgenomics.com/hc/en-us/articles/360022448251-How-to-filter-the-BAM-file-produced-by-10x-pipelines-with-a-list-of-barcodes-
process RENAME_READS {
    tag { "grep CB reads in ${sample_id}" }
    // label "process_low"
    publishDir "${params.outdir}/process_cr/", mode: 'symlink'

    input:
    // sam containing mapped reads
    tuple val(sample_id), path(sam)
    // input bam (output of cellranger)
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.mapped_renamed.bam")

    script:
    """
    # Save the header lines
samtools view -H $sam > SAM_header.sam

# get lines from input bam (containing cb and umi info) that are mapped to ref
awk -F"\t" 'FNR==NR {lines[\$1]; next} \$1 in lines' <(samtools view $sam) <(samtools view $bam) > ${sample_id}_rename_helper.bam
# sort to be able to past later
samtools sort -n -@ 8 -O sam ${sample_id}_rename_helper.bam > ${sample_id}_rename_helper_sorted.sam
samtools sort -n -@ 8 -O sam $sam > ${sample_id}.mapped.sorted.sam
# create header that contains cb and umi, needed for umi_tools count
sed 's/^\\([^\t]*\\).*CB:Z:\\([^\t]*\\)-1.*UB:Z:\\([^\t]*\\).*/\\1\t\\1_\\2_\\3/g' <(samtools view ${sample_id}_rename_helper_sorted.sam) > read_names.txt
# paste new header and rest of sam together
# paste -d '\t' read_names.txt <(awk -v OFS="\t" '{\$1=""; print substr(\$0,2)}' <(samtools view ${sample_id}.mapped.sorted.sam)) > filtered_SAM_body_renamed.sam
# join based on first column, then print all but first column
join -j1 -t "\$(echo -e "\t")" read_names.txt <(samtools view ${sample_id}.mapped.sorted.sam) | awk -v OFS="\t" '{\$1=""; print substr(\$0,2)}' > filtered_SAM_body_renamed.sam

# Combine header and body
cat SAM_header.sam filtered_SAM_body_renamed.sam > ${sample_id}.mapped_renamed.sam
samtools view -bS ${sample_id}.mapped_renamed.sam > ${sample_id}.mapped_renamed.bam
    """
}
