process BENCHMARK_PEAKS {
    tag "$meta.id"
    label 'process_medium'
    label "error_ignore"

    //conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://giorgiagandolfi97/process_1.3:05' :
        'docker.io/giorgiagandolfi97/process_1.3:05' }"

    publishDir { "${params.outdir}/process_simulation/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(peaks), path(cell_rds_files)

    
    output:
    tuple val(meta), path("*.rds"), emit: peak_benchmark_rds
    
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    # The single peaks file for this sample
    peaks <- read.table("${peaks}", header = FALSE, sep = "\\t")

    # All per-cell RDS files (grouped by sample_id, staged as a list)
    #cell_files <- list.files(path = ".", pattern = "fragments_cell_id_.*", full.names = TRUE)
    x <- "${cell_rds_files}"
    files <- strsplit(x, " ")[[1]]
    peak_accessibility_list <- lapply(files, readRDS)
    #cell_data <- lapply(cell_files, readRDS)

    # --- Your analysis logic here ---

    saveRDS(cell_data, file = "${meta.sample_id}_aggregated.rds")
    """

    //template "main_script.R"
}