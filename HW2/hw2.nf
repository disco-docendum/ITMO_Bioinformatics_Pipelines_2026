params.sra_id    = null
params.reads     = null
params.reference = null
params.outdir    = "results"


// 1. Fetch SRA if provided
process FETCH_SRA {
    tag "$sra_id"
    conda 'bioconda::sra-tools=3.0.3'
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
    val sra_id

    output:
    tuple val(sra_id), path("${sra_id}_{1,2}.fastq")

    script:
    """
    fasterq-dump --threads ${task.cpus} --split-files $sra_id
    """
}

// 2. Raw QC
process FASTQC_RAW {
    tag "$sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} -q ${reads[0]} ${reads[1]}
    """
}

// 3. Trim reads
process TRIMMOMATIC {
    tag "$sample_id"
    conda 'bioconda::trimmomatic=0.39'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R?_p.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} $reads \\
        ${sample_id}_R1_p.fastq ${sample_id}_R1_u.fastq \\
        ${sample_id}_R2_p.fastq ${sample_id}_R2_u.fastq \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

// 4. Trimmed QC
process FASTQC_TRIMMED {
    tag "$sample_id"
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} -q ${reads[0]} ${reads[1]}
    """
}

// 5. De novo assembly (if no reference is provided)
process SPADES {
    tag "$sample_id"
    conda 'bioconda::spades=3.15.5 conda-forge::python=3.10'
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    spades.py -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o spades_out
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}

// 6. Map reads
process MAP_READS {
    tag "$sample_id"
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'
    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa index $reference
    bwa mem -t ${task.cpus} $reference ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -bS - | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    """
}

// 7. Plots
process PLOT_COVERAGE {
    tag "$sample_id"
    conda 'bioconda::samtools=1.19 conda-forge::python=3.10 conda-forge::matplotlib=3.8.0'
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_depth.txt")

    script:
    """
    samtools depth $bam > ${sample_id}_depth.txt

    cat <<EOF > plot_coverage.py
    import matplotlib.pyplot as plt

    depths = []
    with open("${sample_id}_depth.txt") as f:
        for line in f:
            # samtools depth format: chr  pos  depth
            depths.append(int(line.strip().split()[2]))

    plt.figure(figsize=(12, 4))
    plt.plot(depths, color='darkred', linewidth=0.5)
    plt.title('Coverage Depth Across Reference: ${sample_id}')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('${sample_id}_coverage.png', dpi=300)
    EOF

    python3 plot_coverage.py
    """
    }

// Main workflow
workflow {
    
    // 1. Initialize
    if (params.sra_id) {
        raw_reads_ch = FETCH_SRA(Channel.of(params.sra_id))
    } else if (params.reads) {
        raw_reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    } else {
        error "Please provide either '--sra_id' or '--reads \"data/*_{1,2}.fastq\"'"
    }

    // 2 & 3 & 4. QC and trimming
    FASTQC_RAW(raw_reads_ch)
    trimmed_reads_ch = TRIMMOMATIC(raw_reads_ch)
    FASTQC_TRIMMED(trimmed_reads_ch)

    // 5. Reference assembly
    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
        mapping_in_ch = trimmed_reads_ch.combine(ref_ch)
    } else {
        assembly_ch = SPADES(trimmed_reads_ch)
        mapping_in_ch = trimmed_reads_ch.join(assembly_ch)
    }

    // 6 & 7. Mapping and plotting
    bam_ch = MAP_READS(mapping_in_ch)
    PLOT_COVERAGE(bam_ch)
}