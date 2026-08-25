process ART_MODERN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::art_modern=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/art_modern:1.3.2--hb45bfb9_0' :
        'biocontainers/art_modern:1.3.2--hb45bfb9_0' }"
        
    publishDir { "${params.outdir}/fastq/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.pe.fastq"), emit: fastq
    path "versions.yml"                         , emit: versions

    
    
    script:
    def fcov     = task.ext.args?.contains('--i-fcov') ? '' : "--i-fcov ${params.fold_coverage}"
    def read_len = task.ext.args?.contains('--read_len') ? '' : "--read_len ${params.read_length}"
    def min_qual = task.ext.args?.contains('--min_qual') ? '' : "--min_qual ${params.min_qual}"
    def args     = task.ext.args ?: ''
    """
    art_modern \\
        --mode template \\
        --lc pe \\
        --i-file ${fasta} \\
        --o-fastq ${meta.id}.pe.fastq \\
        ${fcov} \\
        ${read_len} \\
        ${min_qual} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        art_modern: \$(art_modern --version 2>&1 | head -n1 | sed -e 's/.*version //')
    END_VERSIONS
    """
}
