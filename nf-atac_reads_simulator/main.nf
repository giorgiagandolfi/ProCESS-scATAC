/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Simulate scATAC-seq fragments pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three-step workflow:
    1. Extract FASTA sequences from per-cell BED intervals (bedtools getfasta)
    2. Simulate paired-end reads from those sequences (art_modern)
    3. Split interleaved FASTQ into R1/R2 (bbmap reformat.sh)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ----------------------------
// Process: BEDTOOLS_GETFASTA
// ----------------------------
process BEDTOOLS_GETFASTA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"
        
    publishDir { "${params.outdir}/fasta/${meta.id}" }, mode: 'copy'
    
    
    
    
    input:
    tuple val(meta), path(bed)
    path(fasta)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: fasta
    path "versions.yml"                   , emit: versions
    
    

    script:
    """
    bedtools getfasta -name \\
        -fi ${fasta} \\
        -bed ${bed} \\
        -fo ${meta.id}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//')
    END_VERSIONS
    """
}

// ----------------------------
// Process: ART_MODERN
// ----------------------------
process ART_MODERN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::art_modern=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/art_modern:1.3.2--h4b8a817_0' :
        'biocontainers/art_modern:1.3.2--h4b8a817_0' }"
        
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

// ----------------------------
// Process: BBMAP_REFORMAT
// ----------------------------
process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bbmap=39.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap%3A39.81--h9b5c0a0_1' :
        'biocontainers/bbmap:39.81--he5f23fd_1' }"
    
    publishDir { "${params.outdir}/fastq/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads
    path "versions.yml"                                                             , emit: versions
    
    
    
    script:
    """
    reformat.sh \\
        in=${fastq} \\
        out1=${meta.id}_R1.fastq.gz \\
        out2=${meta.id}_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh 2>&1 | head -n1)
    END_VERSIONS
    """
}



process ALIGN_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/46/46865098b2b8fac5aa4c4200406bb739bf7d06cd9ee9b502b512b9a7f97cbd14/data' :
        'community.wave.seqera.io/library/bowtie2_samtools:3daa9c846bba9a07' }"

    publishDir { "${params.outdir}/aligned_bam/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path(bowtie2_index)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: cb_bam
    path "versions.yml"                    , emit: versions
    

    script:
    def bowtie2_prefix = bowtie2_index[0].name.toString().replaceAll(/\.(rev\.)?\d+\.bt2$/, '')
    
    """
    bowtie2 \\
        --very-sensitive \\
        -X 2000 \\
        -x ${bowtie2_prefix} \\
        -1 ${fastq_1} \\
        -2 ${fastq_2} \\
        --rg-id ${meta.id} \\
        --rg "SM:${meta.id}" \\
    | awk -v CB="${meta.id}" 'BEGIN{OFS="\\t"} /^@/{print; next} {print \$0, "CB:Z:" CB}' \\
    | samtools sort -o ${meta.id}.bam -
    
    samtools index ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n1 | sed 's/.*version //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process SAMTOOLS_MERGE {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    publishDir { "${params.outdir}/merged_bam/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.sample_id}.merged.bam"), path("${meta.sample_id}.merged.bam.bai"), emit: merged_bam
    path "versions.yml"                                                                             , emit: versions

    script:
    """
    samtools merge ${meta.sample_id}.merged.bam --threads ${task.cpus - 1} ${bams}
    samtools index ${meta.sample_id}.merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}



process SINTO_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::sinto=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sinto:0.10.1--pyhdfd78af_0' :
        'biocontainers/sinto:0.10.1--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample_id}.bed"), emit: frag_bed
    path "versions.yml"                           , emit: versions

    script:
    """
    sinto fragments \\
        --bam ${bam} \\
        -f ${meta.sample_id}.bed \\
        -t CB \\
        --use_chrom '.*'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version 2>&1 | head -n1)
    END_VERSIONS
    """
}


process SORT_FRAGMENTS {
    tag "$meta.sample_id"
    label 'process_medium'

    conda "bioconda::sinto=0.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabixpp%3A1.1.0--hd2e4403_0' :
        'biocontainers/sinto:0.10.1--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/fragment_calling/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.sample_id}.fragments.tsv.gz"), emit: frag_gz
    path "versions.yml"                           , emit: versions

    script:
    """
    sort -k1,1 -k2,2n ${bed} > ${meta.sample_id}.fragments.tsv

    bgzip -@ 4 ${meta.sample_id}.fragments.tsv

    tabix -p bed ${meta.sample_id}.fragments.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version 2>&1 | head -n1)
    END_VERSIONS
    """
}

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


// ----------------------------
// Workflow
// ----------------------------
workflow {

    ch_beds = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.cell_id, sample_id: row.sample_id]
            def bed  = file(row.cell_bed, checkIfExists: true)
            [meta, bed]
        }

    ch_fasta = channel.fromPath(params.fasta, checkIfExists: true).first()
    ch_bt2_index = channel.fromPath("${params.bowtie2_index}*.bt2").collect()

    // Step 1: Extract FASTA sequences from BED intervals
    BEDTOOLS_GETFASTA(ch_beds, ch_fasta)

    // Step 2: Simulate paired-end reads from extracted sequences
    ART_MODERN(BEDTOOLS_GETFASTA.out.fasta)

    // Step 3: Split interleaved FASTQ into R1 and R2
    BBMAP_REFORMAT(ART_MODERN.out.fastq)

    // Step 4: Align reads
    ALIGN_BAM(BBMAP_REFORMAT.out.reads, ch_bt2_index)
    
    // Step 5: Group by sample and merge
    ch_bams_grouped = ALIGN_BAM.out.cb_bam
      .map { meta, bam, bai ->
          def group_meta = [sample_id: meta.sample_id]
          [group_meta, bam, bai]
      }
      .groupTuple()

    SAMTOOLS_MERGE(ch_bams_grouped)
    
    // call fragments
    SINTO_FRAGMENTS(SAMTOOLS_MERGE.out.merged_bam)
    SORT_FRAGMENTS(SINTO_FRAGMENTS.out.frag_bed)
    
    CALL_PEAKS(SAMTOOLS_MERGE.out.merged_bam)
}