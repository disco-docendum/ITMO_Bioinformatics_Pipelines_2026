nextflow.enable.dsl=2

// Import processes and workflows from modules file and nf-core
include { FETCH_SRA; FASTQC; TRIMMOMATIC; TRIMMED_QC_WF; FETCH_REFERENCE; SPADES; MAP_READS; PLOT_STATS } from './modules.nf'

include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL } from './modules/nf-core/bcftools/call/main'

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

    // 4. Mapping
    bam_ref_ch = MAP_READS(mapping_in_ch)

    // 5. Variant calling (The nf-core way)
    ch_split = bam_ref_ch.multiMap { sample_id, bam, reference -> 
        bam_input:   tuple( [id: sample_id, single_end: false], bam, [], [] ) 
        fasta_input: tuple( [id: 'reference'], reference, [] )
        plot_input:  tuple( sample_id, bam )
    }

    // Run mpileup
    mpileup_out = BCFTOOLS_MPILEUP(ch_split.bam_input, ch_split.fasta_input.first(), false)

    // 6. Coverage and variant plotting
    clean_vcf_ch = mpileup_out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }
    
    plot_in_ch = ch_split.plot_input.join(clean_vcf_ch)
    PLOT_STATS(plot_in_ch)
}
