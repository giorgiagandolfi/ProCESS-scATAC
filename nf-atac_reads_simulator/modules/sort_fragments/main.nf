process SORT_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::sinto=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabixpp%3A1.1.0--hd2e4403_0' :
        'biocontainers/sinto:0.10.1--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.sample_id}.sorted.bed.gz"), emit: frag_gz
    path "versions.yml"                           , emit: versions

    script:
    """
    sort -k1,1 -k2,2n ${bed} > ${meta.sample_id}.sorted.bed

    bgzip -@ 4 ${meta.sample_id}.sorted.bed

    tabix -p bed ${meta.sample_id}.sorted.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version 2>&1 | head -n1)
    END_VERSIONS
    """
}