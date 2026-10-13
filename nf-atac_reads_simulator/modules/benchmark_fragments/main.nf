process BENCHMARK_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_medium'
    label "error_ignore"

    //conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://giorgiagandolfi97/my_genomics:v0' :
        'docker.io/giorgiagandolfi97/my_genomics:v0' }"

    publishDir { "${params.outdir}/benchmarking/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fragments), path(cell_peak_txt_files)
    

    
    output:
    tuple val(meta), path("*.rds"), emit: fragment_benchmark_rds
    tuple val(meta), path("*.pdf"), emit: fragment_benchmark_report_pdf
    
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    template "main_script.R"
}