process CALL_PEAKS_HMMRATAC {
    tag "$meta.sample_id"
    label 'process_high'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs3:3.0.4--py313hda738de_0' :
        'quay.io/biocontainers/macs3:3.0.4--py313hda738de_0' }"

    publishDir {
        "${params.outdir}/call_peaks/hmmratac/${meta.sample_id}"
    }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peak
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    macs3 hmmratac \\
        --input ${bam} \\
        --format BAMPE \\
        --name ${prefix} \\
        ${args}

    # Extract genomic coordinates for benchmarking
    # from the HMMRATAC accessible-region output.
    if [ -f "${prefix}_accessible_regions.narrowPeak" ]; then
        mv "${prefix}_accessible_regions.narrowPeak" \\
           "${prefix}.hmmratac.narrowPeak"
    else
        echo "ERROR: HMMRATAC accessible-region output not found" >&2
        exit 1
    fi

    cut -f1-3 "${prefix}.hmmratac.narrowPeak" \\
        > "${prefix}.hmmratac.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version 2>&1 | sed 's/macs3 //')
    END_VERSIONS
    """
}
