process SORT_BAM_BY_NAME {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"


    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.namesorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    samtools sort \\
        -n \\
        -@ ${task.cpus} \\
        -o ${meta.sample_id}.namesorted.bam \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}