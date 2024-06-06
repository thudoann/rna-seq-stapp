nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */
params.single_end_reads = "Data/${params.study_id}/*_pass.fastq"
params.transcriptome = "t_indxs/pao1_cnda.fa"
params.outdir = "Output/${params.study_id}"
params.kmer_size = 15  // Default value, can be overridden by user input
params.study_id = "default_study"  // Study identifier, should be provided by user input

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         single-end reads: ${params.single_end_reads}
         outdir       : ${params.outdir}
         kmer size    : ${params.kmer_size}
         """
         .stripIndent()

/*
 * define the `INDEX` process that creates a binary index
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
 * Define TrimGalore process for trimming single-end reads
 */
process TRIMGALORE_SINGLE {
    cpus 2
    tag "TrimGalore on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq")

    script:
    """
    trim_galore -o . ${reads}
    if [ "${reads.baseName}_trimmed.fq" != "${sample_id}_trimmed.fq" ]; then
        mv ${reads.baseName}_trimmed.fq ${sample_id}_trimmed.fq
    fi
    """
}

/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process QUANT_SINGLE {
    cpus 4
    tag "Quantification on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/quant", mode: 'copy'

    input:
    each index
    tuple val(sample_id), path(reads)

    output:
    path(sample_id)

    script:
    """
    salmon quant --threads ${task.cpus} -i ${index} --libType A -r ${reads} --validateMappings -o ${sample_id}
    """
}

/*
 * Run fastQC to check quality of reads files
 */
process FASTQC_SINGLE {
    cpus 2
    tag "FASTQC on $sample_id (single-end)"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
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
    single_end_reads_ch = Channel.fromPath("Data/${params.study_id}/fastq/*_pass.fastq", checkIfExists: true)
                            .map { file -> 
                                println "Found file: ${file}"
                                tuple(file.simpleName.replace('_pass', ''), file) 
                            }

    transcriptome_ch = Channel.fromPath(params.transcriptome, checkIfExists: true)

    index_ch = INDEX(transcriptome_ch)

    trimmed_single_end_reads_ch = TRIMGALORE_SINGLE(single_end_reads_ch)

    fastqc_single_ch = FASTQC_SINGLE(trimmed_single_end_reads_ch)
    quant_single_ch = QUANT_SINGLE(index_ch, trimmed_single_end_reads_ch)

    multiqc_ch = MULTIQC(quant_single_ch.mix(fastqc_single_ch).collect())
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone! Open the following report in your browser --> ${params.outdir}/multiqc_report.html\n" : "Oops .. something went wrong")
}
