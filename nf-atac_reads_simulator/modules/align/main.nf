process ALIGN_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46865098b2b8fac5aa4c4200406bb739bf7d06cd9ee9b502b512b9a7f97cbd14/data' :
        'community.wave.seqera.io/library/bowtie2_samtools:3daa9c846bba9a07' }"

    publishDir { "${params.outdir}/aligned_bam/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path(bowtie2_index)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: cb_bam
    path "versions.yml"                    , emit: versions
    

    script:
    def bowtie2_prefix = bowtie2_index[0].name.toString().replaceAll(/\.(rev\.)?\d+\.bt2$/, '')
    
    """
    bowtie2 \\
        --very-sensitive \\
        -X 2000 \\
        -x ${bowtie2_prefix} \\
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        --rg-id ${meta.id} \\
        --rg "SM:${meta.id}" \\
    | awk -v CB="${meta.id}" 'BEGIN{OFS="\\t"} /^@/{print; next} {print \$0, "CB:Z:" CB}' \\
    | samtools sort -o ${meta.id}.bam -
    
    samtools index ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n1 | sed 's/.*version //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}