process CALL_PEAKS {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs3:3.0.4--py313hda738de_0' :
        'biocontainers/macs3:3.0.4--py313hda738de_0' }" //check this later on

    publishDir { "${params.outdir}/call_peaks/macs3/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*.xls")                   , emit: xls
    path "versions.yml"                  , emit: versions
    tuple val(meta), path("*.gappedPeak"), optional:true, emit: gapped
    tuple val(meta), path("*.bed")       , optional:true, emit: bed
    tuple val(meta), path("*.bdg")       , optional:true, emit: bdg

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    //def format    = meta.single_end ? 'BAM' : 'BAMPE'
    //def control   = controlbam ? "--control $controlbam" : ''
    """
    macs3 \\
        callpeak \\
        ${args} \\
        -g 2.7e9 \\
        --format 'BAMPE' \\
        --name $prefix \\
        --treatment $bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version 2>&1 | sed 's/macs3 //')
    END_VERSIONS
    """
}