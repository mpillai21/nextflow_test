#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// set param values
params.reads = "$baseDir/fastq/*_{R1,R2}.fastq.gz"
params.outdir = "$baseDir/results"

// Process for trimming reads using Trimmomatic
process TrimReads {
    tag "${pair_id}"

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(pair_id), path(read1), path(read2)

    output:
    tuple val(pair_id), path("${pair_id}_R1.trimmed.fastq.gz"), path("${pair_id}_R2.trimmed.fastq.gz")

    script:
    """
    trimmomatic PE -phred33 \\
        $read1 $read2 \\
        ${pair_id}_R1.trimmed.fastq.gz ${pair_id}_R1.unpaired.fastq.gz \\
        ${pair_id}_R2.trimmed.fastq.gz ${pair_id}_R2.unpaired.fastq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process for assembling trimmed reads using SPAdes
process AssembleFasta {
    tag "${pair_id}"
    
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    tuple val(pair_id), path(r1_trimmed), path(r2_trimmed)

    output:
    path("${pair_id}.assembly")

    script:
    """
    spades.py -1 $r1_trimmed -2 $r2_trimmed -o ${pair_id}.assembly
    """
}

// Workflow block
workflow {
    // Create Channel from paired-end read files
    read_pairs = Channel.fromFilePairs(params.reads, size: 2, flat: true)

    // Run TrimReads process
    trimmed_reads_ch = read_pairs
        | TrimReads

    // Run AssembleFasta process
    assembled_genome_ch = trimmed_reads_ch
        | AssembleFasta

    // Display completion message for each assembled fasta file
    assembled_genome_ch.view { "Assembly completed: ${it}" }
}

