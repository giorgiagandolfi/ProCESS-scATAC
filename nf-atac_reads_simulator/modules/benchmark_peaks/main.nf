process BENCHMARK_PEAKS {
    tag "$meta.sample_id"
    label 'process_medium'
    label "error_ignore"

    //conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://giorgiagandolfi97/my_genomics:v0' :
        'docker.io/giorgiagandolfi97/my_genomics:v0' }"

    publishDir { "${params.outdir}/benchmarking/peak_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(peaks), path(cell_peak_txt_files), path(cell_frags_rds_files)

    
    output:
    tuple val(meta), path("*.rds"), emit: peak_benchmark_rds
    
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template "main_script.R"
}
