process PROCESS_SIMULATE_OUT_PEAK_FRAGMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/macs3:3.0.4--py313hda738de_0' :
        'biocontainers/macs3:3.0.4--py313hda738de_0' }" //check this later on

    publishDir { "${params.outdir}/process_simulation/out_peak/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(in_peak_frg)
    path(gap_file)
    path(centromere_file)
    path(frag_dist_out)
    
    output:
    tuple val(meta), path("*.rds"), emit: out_peak_frg
    tuple val(meta), path("*.fasta")                   , emit: out_peak_fasta
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    template "main_script.R"
}
