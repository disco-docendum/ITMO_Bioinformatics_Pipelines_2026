// 1. Fetch SRA
process FETCH_SRA {
    tag "$sra_id"
    conda 'bioconda::sra-tools=3.0.3'
    container 'quay.io/biocontainers/sra-tools:3.0.3--h860ed49_0'
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input: val sra_id
    output: tuple val(sra_id), path("${sra_id}_{1,2}.fastq")

    script: "fasterq-dump --threads ${task.cpus} --split-files $sra_id"
}

// 2. QC
process FASTQC {
    tag "$sample_id"
    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc_${reads_type}", mode: 'copy'

    input:
    val reads_type
    tuple val(sample_id), path(reads)

    output: path "*_fastqc.{zip,html}"

    script: "fastqc -t ${task.cpus} -q ${reads[0]} ${reads[1]}"
}

// 3. Trimming
process TRIMMOMATIC {
    tag "$sample_id"
    conda 'bioconda::trimmomatic=0.39'
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: "*_p.fastq"

    input: tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("${sample_id}_R?_p.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} $reads \\
        ${sample_id}_R1_p.fastq ${sample_id}_R1_u.fastq \\
        ${sample_id}_R2_p.fastq ${sample_id}_R2_u.fastq \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

workflow TRIMMED_QC_WF {
    take: reads
    main: FASTQC('trimmed', reads)
    emit: FASTQC.out
}

// 4. NEW: Fetch reference genome from NCBI!
process FETCH_REFERENCE {
    tag "$ref_id"
    conda 'bioconda::entrez-direct=16.2'
    container 'quay.io/biocontainers/entrez-direct:16.2--he881be0_1'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val ref_id

    output:
    path "${ref_id}.fasta"

    script:
    """
    efetch -db nuccore -id $ref_id -format fasta > ${ref_id}.fasta
    """
}


// 5. De novo assembly
process SPADES {
    tag "$sample_id"
    conda 'bioconda::spades=3.15.5 conda-forge::python=3.10'
    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input: tuple val(sample_id), path(reads)
    output: tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    spades.py -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o spades_out
    cp spades_out/contigs.fasta ${sample_id}_contigs.fasta
    """
}

// 6. Mapping
process MAP_READS {
    tag "$sample_id"
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.19'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79fd4d6:a208d0e2e92c2a4ddba2417b1fb9e19c3b8cb4e6-0'
    publishDir "${params.outdir}/mapped", mode: 'copy'

    input: tuple val(sample_id), path(reads), path(reference)
    output: tuple val(sample_id), path("${sample_id}.sorted.bam"), path(reference)

    script:
    """
    bwa index $reference
    bwa mem -t ${task.cpus} $reference ${reads[0]} ${reads[1]} | samtools view -@ ${task.cpus} -bS - | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    """
}

// 8. Plots (now both for coverage and variants!)
process PLOT_STATS {
    tag "$sample_id"
    conda 'bioconda::samtools=1.19 conda-forge::python=3.10 conda-forge::matplotlib=3.8.0'
    publishDir "${params.outdir}/plots", mode: 'copy'

    input: 
    tuple val(sample_id), path(bam), path(vcf)

    output: 
    tuple val(sample_id), path("${sample_id}_coverage.png"), path("${sample_id}_variants.png"), path("${sample_id}_depth.txt")

script:
    """
    samtools depth ${bam} > ${sample_id}_depth.txt

    cat <<EOF > plot_data.py
    import matplotlib.pyplot as plt
    import gzip

    # 1. Plot coverage
    depths = []
    with open("${sample_id}_depth.txt") as f:
        for line in f:
            depths.append(int(line.strip().split()[2]))

    plt.figure(figsize=(12, 4))
    plt.plot(depths, color='darkred', linewidth=0.5)
    plt.title('Coverage depth across reference: ${sample_id}')
    plt.xlabel('Position')
    plt.ylabel('Depth')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('${sample_id}_coverage.png', dpi=300)
    plt.close()

    # 2. Plot variants coordinates
    variant_positions = []
    
    with gzip.open("${vcf}", "rt") as f:
        for line in f:
            # Skip VCF header lines
            if not line.startswith('#'):
                parts = line.split('\t')
                variant_positions.append(int(parts[1]))

    plt.figure(figsize=(12, 2))
    if variant_positions:
        plt.eventplot(variant_positions, color='black', linewidths=1.5)

    plt.title('Variant locations across reference: ${sample_id}')
    plt.xlabel('Genomic position')
    plt.yticks([])

    if depths:
        plt.xlim(0, len(depths))

    plt.tight_layout()
    plt.savefig('${sample_id}_variants.png', dpi=300)
    EOF

    python3 plot_data.py
    """
}
