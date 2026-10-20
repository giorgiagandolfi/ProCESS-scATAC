process BENCHMARK_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_high_memory'
    label "error_ignore"

    //conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://giorgiagandolfi97/my_genomics:v0' :
        'docker.io/giorgiagandolfi97/my_genomics:v0' }"

    publishDir { "${params.outdir}/benchmarking/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fragments), path(cell_frags_rds_files)
    

    
    output:
    tuple val(meta), path("*_fragments_metrics.rds"), emit: fragment_benchmark_rds
    tuple val(meta), path("*_fragments_metrics_plots.rds"), emit: fragment_benchmark_report_rds
    
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template "main_script.R"
}
