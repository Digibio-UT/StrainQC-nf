#!/usr/bin/env nextflow
/*
 * StrainQC-nf — entry point.
 *
 * Usage:
 *   nextflow run . -profile docker \
 *       --input samplesheet.csv --reference ref.gbk \
 *       --variant_caller gatk --ploidy 1
 *
 * Samplesheet CSV columns: sample,fastq_1,fastq_2   (fastq_2 optional for single-end)
 */
nextflow.enable.dsl = 2

include { STRAINQC } from './workflows/strainqc.nf'

workflow {
    if (!params.input)     error "Provide --input <samplesheet.csv>"
    if (!params.reference) error "Provide --reference <genome.gbk|.fasta>"

    // Parse samplesheet → [ meta, [reads] ]
    ch_reads = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta  = [ id: row.sample ]
            def reads = row.fastq_2 ? [ file(row.fastq_1, checkIfExists: true),
                                        file(row.fastq_2, checkIfExists: true) ]
                                    : [ file(row.fastq_1, checkIfExists: true) ]
            [ meta, reads ]
        }

    ch_reference = file(params.reference, checkIfExists: true)

    STRAINQC(ch_reads, ch_reference)
}

workflow.onComplete {
    log.info ( workflow.success
        ? "✅ StrainQC-nf finished. Results in ${params.outdir}/"
        : "❌ StrainQC-nf failed — see .nextflow.log" )
}
