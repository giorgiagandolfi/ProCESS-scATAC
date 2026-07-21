process PROCESS_SIMULATE_FRAGMENT {
    tag "$meta.id"
    label 'process_medium'
    label "error_ignore"

    //conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://giorgiagandolfi97/process_1.3:05' :
        'docker.io/giorgiagandolfi97/process_1.3:05' }"

    publishDir { "${params.outdir}/process_simulation/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(sample_forest), path(pyhlo_forest)
    path(activity_list_file)
    path(peak_list_file)
    path(fragment_size_distribution)
    path(encode_blacklist)
    path(gap_file)
    path(centromere_file)
    path(fragm_len_dist_out)
    
    output:
    tuple val(meta), path("*.rds"), emit: all_fragments_rds
    tuple val(meta), path("*.fasta")                   , emit: all_fragments_fasta
    tuple val(meta), path("*.txt")                   , emit: peak_fragment_mapping
    
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    template "main_script.R"
}