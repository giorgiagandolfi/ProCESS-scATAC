process CALL_PEAKS_GENRICH {
    tag "$meta.sample_id"
    label 'process_high'

    conda "bioconda::genrich=0.6.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genrich:0.6.2--h118bc1c_0' :
        'quay.io/biocontainers/genrich:0.6.2--h118bc1c_0' }"

    publishDir { "${params.outdir}/call_peaks/genrich/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peak
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    Genrich \\
        ${args} \\
        -t ${bam} \\
        -o ${prefix}.genrich.narrowPeak \\
        -j \\
        -q 0.05 \\
        -v

    cut -f1-3 ${prefix}.genrich.narrowPeak \\
        > ${prefix}.genrich.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genrich: \$(Genrich -V 2>&1 | head -n 1)
    END_VERSIONS
    """
}