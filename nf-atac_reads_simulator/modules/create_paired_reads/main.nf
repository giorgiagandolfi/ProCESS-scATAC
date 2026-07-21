process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bbmap=39.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap%3A39.81--h9b5c0a0_1' :
        'biocontainers/bbmap:39.81--he5f23fd_1' }"
    
    publishDir { "${params.outdir}/fastq/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads
    path "versions.yml"                                                             , emit: versions
    
    
    
    script:
    """
    reformat.sh \\
        in=${fastq} \\
        out1=${meta.id}_R1.fastq.gz \\
        out2=${meta.id}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh 2>&1 | head -n1)
    END_VERSIONS
    """
}