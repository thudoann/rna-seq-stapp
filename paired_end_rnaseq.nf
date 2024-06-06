nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */
params.paired_end_reads = "Data/${params.study_id}/fastq/*_pass_{1,2}.fastq"
params.transcriptome = "t_indxs/pao1_cnda.fa"
params.outdir = "Output/${params.study_id}"
params.kmer_size = 15  // Default value, can be overridden by user input
params.study_id = "default_study"  // Study identifier, should be provided by user input

log.info """
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome     : ${params.transcriptome}
         paired-end reads  : ${params.paired_end_reads}
         outdir            : ${params.outdir}
         kmer size         : ${params.kmer_size}
         """
         .stripIndent()

/*
 * Define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
    cpus 4
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads ${task.cpus} -t ${transcriptome} -i index -k ${params.kmer_size}
    """
}

/*
 * Define TrimGalore process for trimming paired-end reads
 */
process TRIMGALORE_PAIRED {
    cpus 2
    memory '8 GB'
    tag "TrimGalore on $pair_id (paired-end)"
    publishDir "${params.outdir}/${pair_id}/trimmed", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("${reads[0].baseName}_trimmed.fq"), path("${reads[1].baseName}_trimmed.fq")

    script:
    """
    set -e
    echo "Running trim_galore on ${reads[0]} and ${reads[1]}"
    trim_galore --paired -o . ${reads[0]} ${reads[1]}
    if [ ! -f ${reads[0].baseName}_val_1.fq ] || [ ! -f ${reads[1].baseName}_val_2.fq ]; then
        echo "Trim Galore output files not found."
        exit 1
    fi
    mv ${reads[0].baseName}_val_1.fq ${reads[0].baseName}_trimmed.fq
    mv ${reads[1].baseName}_val_2.fq ${reads[1].baseName}_trimmed.fq
    """
}


/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process QUANT_PAIRED {
    cpus 4
    tag "Quantification on $pair_id (paired-end)"
    publishDir "${params.outdir}/${pair_id}/quant", mode: 'copy'

    input:
    each index
    tuple val(pair_id), path(read1), path(read2)

    output:
    path(pair_id)

    script:
    """
    salmon quant --threads ${task.cpus} -i ${index} --libType A -1 ${read1} -2 ${read2} --validateMappings -o ${pair_id}
    """
}

/*
 * Run FastQC to check the quality of read files
 */
process FASTQC_PAIRED {
    cpus 2
    tag "FASTQC on $pair_id (paired-end)"
    publishDir "${params.outdir}/${pair_id}/fastqc", mode: 'copy'

    input:
    tuple val(pair_id), path(read1), path(read2)

    output:
    path("fastqc_${pair_id}_logs")

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${read1} ${read2}
    """
}

/*
 * Create a report using multiQC for the quantification
 * and fastqc processes
 */
process MULTIQC {
    cpus 2
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path('*')

    output:
    path('multiqc_report.html')

    script:
    """
    multiqc .
    """
}


workflow {
    paired_end_reads_ch = Channel.fromFilePairs(params.paired_end_reads, checkIfExists: true)

    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)

    index_ch = INDEX(transcriptome_ch)

    trimmed_paired_end_reads_ch = TRIMGALORE_PAIRED(paired_end_reads_ch)

    fastqc_paired_ch = FASTQC_PAIRED(trimmed_paired_end_reads_ch)
    quant_paired_ch = QUANT_PAIRED(index_ch, trimmed_paired_end_reads_ch)

    multiqc_ch = MULTIQC(quant_paired_ch.mix(fastqc_paired_ch).collect())
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone! Open the following report in your browser --> ${params.outdir}/multiqc_report.html\n" : "Oops .. something went wrong")
}
