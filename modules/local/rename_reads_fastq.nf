// add cell barcode and UMI to read name
// necessary for BAM from cellranger or STARsolo as input
// performed at this point only on aligned sequences because the process is very slow on large files
process RENAME_READS_FASTQ {
    label "process_high_sc"
    tag "$sample_id"

    input:
        // sam containing mapped reads, bam from cellranger
        tuple val(sample_id), path(fastq), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed_renamed.fastq")

    script:
        """
        # get unmapped reads using many cores
        samtools view -@ ${task.cpus} -f 4 -d CB $bam > unmapped.sam

        # extract fastq
        zcat $fastq > ${sample_id}.trimmed.fastq

        # get lines from input bam (containing cb and umi info) that are mapped to ref
        # create header that contains cb and umi, needed for umi_tools count
        # first column is original header to merge on later
        # piping awk output into sed is faster and uses less CPU and IO and only marginally more memory
        awk -F"\t" 'FNR==NR {lines[\$1]; next} \$1 in lines' <(awk 'NR%4==1{print}' ${sample_id}.trimmed.fastq | sed 's/^@//g') <(cat unmapped.sam) |\\
                parallel -j ${task.cpus} --pipe 'sed "s/^\\([^\t]*\\).*CB:Z:\\([^\t]*\\)-1.*UB:Z:\\([^\t]*\\).*/\\1\t\\1_\\2_\\3/g"' |\\
                sort -S 2G -k1,1 | sed 's/^/@/g' > read_names.txt

        # join 4 rows of a fastq entry into one row
        paste -d "\t"  - - - - < ${sample_id}.trimmed.fastq > ${sample_id}.trimmed.txt

        # join based on first column, then print all but first column
        # split rows into fastq format
        join -j1 -t "\$(echo -e "\t")" read_names.txt <(sort -S 2G -k1,1 ${sample_id}.trimmed.txt) |\\
                awk -v OFS="\t" '{\$1=""; print substr(\$0,2)}' |\\
                sed 's/^/@/g' |\\
                sed 's/\t/\n/g' > ${sample_id}.trimmed_renamed.fastq

        rm ${sample_id}.trimmed.fastq
        rm ${sample_id}.trimmed.txt
        rm unmapped.sam
        """
}
