process CALL_PEAKS_HOMER {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::homer=5.1"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:5.1--pl5321hc52dbad_1' :
        'quay.io/biocontainers/homer:5.1--pl5321hc52dbad_1' }"

    publishDir { "${params.outdir}/call_peaks/homer/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.peaks.txt"), emit: peak
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    makeTagDirectory \\
        ${prefix}_tagdir \\
        ${bam}

    findPeaks \\
        ${prefix}_tagdir \\
        -style dnase \\
        ${args} \\
        -o ${prefix}.homer.peaks.txt

    pos2bed.pl ${prefix}.homer.peaks.txt \\
        > ${prefix}.homer.full.bed

    cut -f1-3 ${prefix}.homer.full.bed \\
        > ${prefix}.homer.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: "5.1"
    END_VERSIONS
    """
}