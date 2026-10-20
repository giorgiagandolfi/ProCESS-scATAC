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
include { CALL_PEAKS_GENRICH }       from "./modules/call_peaks_genrich/main.nf"
include { SORT_BAM_BY_NAME }         from "./modules/sort_by_name/main.nf"  
include { CALL_PEAKS_HOMER }         from "./modules/call_peaks_homer/main.nf"
include { CALL_PEAKS_HMMRATAC }       from "./modules/call_peaks_hmmratac/main.nf"
// GENERAL BENCHMARKING
include { BENCHMARK_PEAKS }               from "./modules/benchmark_peaks/main.nf"
include { BENCHMARK_FRAGMENTS }               from "./modules/benchmark_fragments/main.nf"

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

    // One forest per sample; verify that all cells agree on its path
    ch_sample_forests = ch_cells
        .map { meta, sample_forest, phylo_forest ->
            tuple(meta.sample_id, sample_forest)
        }
        .groupTuple()
        .map { sample_id, forests ->
            def unique_forests = forests.unique()
    
            if (unique_forests.size() != 1) {
                error "Sample ${sample_id} has multiple sample forest paths: ${unique_forests}"
            }
    
            tuple(sample_id, unique_forests[0])
        }
    SAMTOOLS_MERGE(ch_bams_grouped)

    SINTO_FRAGMENTS(SAMTOOLS_MERGE.out.merged_bam)
    SORT_FRAGMENTS(SINTO_FRAGMENTS.out.frag_bed)

    
    CALL_PEAKS(SAMTOOLS_MERGE.out.merged_bam)
    SORT_BAM_BY_NAME(SAMTOOLS_MERGE.out.merged_bam)

    CALL_PEAKS_GENRICH(
        SORT_BAM_BY_NAME.out.bam
    )
    // CALL_PEAKS_HMMRATAC(SAMTOOLS_MERGE.out.merged_bam)
    // CALL_PEAKS_HOMER(SAMTOOLS_MERGE.out.merged_bam)
    
    // Group per-cell RDS files by sample_id
    ch_frags_rds_grouped = PROCESS_SIMULATE_FRAGMENT.out.all_fragments_rds
        .map { meta, rds ->
                def group_meta = [sample_id: meta.sample_id]
                [group_meta, rds]
            }
            .groupTuple()
    ch_peaks_txt_grouped = PROCESS_SIMULATE_FRAGMENT.out.peak_fragment_mapping
        .map { meta, rds ->
                def group_meta = [sample_id: meta.sample_id]
                [group_meta, rds]
            }
            .groupTuple()
    ch_aggregate_input = CALL_PEAKS.out.peak
        .join(ch_peaks_txt_grouped)
        .join(ch_frags_rds_grouped)
        // Use sample_id as the forest join key
        .map { meta, peaks_file, peak_mappings, fragments_rds ->
            tuple(
                meta.sample_id,
                meta,
                peaks_file,
                peak_mappings,
                fragments_rds
            )
        }
        .join(
            ch_sample_forests,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map {
            sample_id, meta, peaks_file, peak_mappings,
            fragments_rds, sample_forest ->
    
            tuple(
                meta,
                peaks_file,
                peak_mappings,
                fragments_rds,
                sample_forest
            )
        }
    

        
    ch_aggregate_input_bp = SORT_FRAGMENTS.out.frag_gz
        .join(ch_frags_rds_grouped)
        
    // Shape: [ [sample_id: X], peaks_file, [rds1, rds2, ...] ]
    ch_aggregate_input.view()
    BENCHMARK_PEAKS(
        ch_aggregate_input,
        ch_accessibility_score,
        ch_peak_file
    )
    BENCHMARK_FRAGMENTS(ch_aggregate_input_bp)
}
