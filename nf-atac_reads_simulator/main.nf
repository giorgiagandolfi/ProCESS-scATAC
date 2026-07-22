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

include { PROCESS_SIMULATE_FRAGMENT } from "./modules/process_simulate_fragments/main.nf"
include { ART_MODERN }               from "./modules/simulate_reads/main.nf"
include { BBMAP_REFORMAT }           from "./modules/create_paired_reads/main.nf"
include { ALIGN_BAM }                from "./modules/align/main.nf"
include { SAMTOOLS_MERGE }           from "./modules/bam_merge/main.nf"
include { SINTO_FRAGMENTS }          from "./modules/call_fragments/main.nf"
include { SORT_FRAGMENTS }           from "./modules/sort_fragments/main.nf"
include { CALL_PEAKS }               from "./modules/call_peaks/main.nf"

workflow {

    // ----------------------------
    // Build the BAM channel depending on entry point
    // ----------------------------
    if (params.step == 'merge_bams') {
        //
        // Start from merge: read a CSV of pre-aligned BAMs grouped by sample
        // CSV columns: sample_id, bam, bai
        //
        
        ch_bams_grouped = channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                // row.cell_id is available but simply not used here
                def meta = [sample_id: row.sample_id]
                [meta, file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true)]
            }
            .groupTuple()

    } else {
        //
        // Full pipeline: simulate fragments ??? reads ??? align
        //
        ch_cells = channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.cell_id, sample_id: row.sample_id]
                def sample_forest_path = file(row.sample_forest_path, checkIfExists: true)
                def phylo_forest_path  = file(row.phylo_forest_path, checkIfExists: true)
                [meta, sample_forest_path, phylo_forest_path]
            }

        ch_fasta             = channel.fromPath(params.fasta, checkIfExists: true).first()
        ch_accessibility_score = channel.fromPath(params.a_score, checkIfExists: true).first()
        ch_peak_file         = channel.fromPath(params.peak_file, checkIfExists: true).first()
        ch_frag_size_dist    = channel.fromPath(params.fragment_size, checkIfExists: true).first()
        ch_blacklist         = channel.fromPath(params.blacklist_file, checkIfExists: true).first()
        ch_gap_file          = channel.fromPath(params.gaps_file, checkIfExists: true).first()
        ch_centromere_file   = channel.fromPath(params.cytoband_file, checkIfExists: true).first()
        ch_fragm_len_dist_out = channel.fromPath(params.fragment_size_out, checkIfExists: true).first()
        ch_bt2_index         = channel.fromPath("${params.bowtie2_index}*.bt2").collect()

        // Step 1: Simulate fragments
        PROCESS_SIMULATE_FRAGMENT(
            ch_cells, ch_accessibility_score, ch_peak_file,
            ch_frag_size_dist, ch_blacklist, ch_gap_file,
            ch_centromere_file, ch_fragm_len_dist_out
        )

        // Step 2: Simulate paired-end reads
        ART_MODERN(PROCESS_SIMULATE_FRAGMENT.out.all_fragments_fasta)

        // Step 3: Split interleaved FASTQ into R1/R2
        BBMAP_REFORMAT(ART_MODERN.out.fastq)

        // Step 4: Align reads
        ALIGN_BAM(BBMAP_REFORMAT.out.reads, ch_bt2_index)

        ALIGN_BAM.out.cb_bam
            .map { meta, bam, bai ->
                def pub_dir = "${params.outdir}/aligned_bam/${meta.id}"
                "${meta.id},${meta.sample_id},${pub_dir}/${bam.name},${pub_dir}/${bai.name}"
            }
            .collectFile(
                name: 'aligned_bams.csv',
                storeDir: "${params.outdir}/csv",
                seed: 'cell_id,sample_id,bam,bai',
                newLine: true,
                sort: true
            )

        // Group by sample for merge
        ch_bams_grouped = ALIGN_BAM.out.cb_bam
            .map { meta, bam, bai ->
                def group_meta = [sample_id: meta.sample_id]
                [group_meta, bam, bai]
            }
            .groupTuple()
    }

    SAMTOOLS_MERGE(ch_bams_grouped)

    SINTO_FRAGMENTS(SAMTOOLS_MERGE.out.merged_bam)
    SORT_FRAGMENTS(SINTO_FRAGMENTS.out.frag_bed)

    CALL_PEAKS(SAMTOOLS_MERGE.out.merged_bam)
}