nextflow.enable.dsl=2

// Import processes and workflows from modules file
include { FETCH_SRA; FASTQC; TRIMMOMATIC; TRIMMED_QC_WF; FETCH_REFERENCE; SPADES; MAP_READS; PLOT_STATS; BCFTOOLS_CALL } from './modules.nf'

workflow {
    
    // 1. Initialize
    if (params.sra_id) {
        raw_reads_ch = FETCH_SRA(Channel.of(params.sra_id))
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Please provide either '--sra_id' or '--reads \"data/*_{1,2}.fastq\"'"
    }

    // 2. QC and trimming
    FASTQC('raw', raw_reads_ch)
    trimmed_reads_ch = TRIMMOMATIC(raw_reads_ch)
    TRIMMED_QC_WF(trimmed_reads_ch)

// 3. Reference assignment, ncbi download, or assembly
    if (params.reference) {
        // Option A: provided local FASTA file
        ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
        mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    } else if (params.ref_id) {
        // Option B: download reference from NCBI dynamically
        ref_ch = FETCH_REFERENCE(params.ref_id)
        mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    } else {
        // Option C: no reference provided, run de novo assembly
        assembly_ch = SPADES(trimmed_reads_ch)
        mapping_in_ch = trimmed_reads_ch.join(assembly_ch)
    }

    // 4. Mapping, output: [sample_id, bam, reference]
    bam_ref_ch = MAP_READS(mapping_in_ch)

    // 5. Variant Calling, output: [sample_id, vcf]
    vcf_ch = BCFTOOLS_CALL(bam_ref_ch)

    // 6. Plotting the results (coverage & variants)
    bam_only_ch = bam_ref_ch.map { sample_id, bam, reference -> tuple(sample_id, bam) }
    
    // Using .join() to combine the BAM and VCF channels based on the sample_id key, output: [sample_id, bam, vcf]
    plot_in_ch = bam_only_ch.join(vcf_ch)
    
    // Pass the combined channel to the plotting process
    PLOT_STATS(plot_in_ch)
}