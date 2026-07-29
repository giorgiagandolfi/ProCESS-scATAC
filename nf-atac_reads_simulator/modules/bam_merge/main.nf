process SAMTOOLS_MERGE {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    publishDir { "${params.outdir}/merged_bam/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.sample_id}.merged.bam"), path("${meta.sample_id}.merged.bam.bai"), emit: merged_bam
    path "versions.yml"                                                                             , emit: versions

    script:
    """
    samtools merge ${meta.sample_id}.merged.bam --threads ${task.cpus - 1} ${bams}
    samtools index ${meta.sample_id}.merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
