process SINTO_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::sinto=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sinto:0.10.1--pyhdfd78af_0' :
        'biocontainers/sinto:0.10.1--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.bed"), emit: frag_bed
    path "versions.yml"                           , emit: versions

    script:
    """
    sinto fragments \\
        --bam ${bam} \\
        -f ${meta.sample_id}.bed \\
        -t CB \\
        --use_chrom '.*'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version 2>&1 | head -n1)
    END_VERSIONS
    """
}