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

    // Parse samplesheet → [ meta, [reads] ]. Read paths may be absolute, a URL, or
    // relative to the samplesheet's own directory (keeps bundled test fixtures portable).
    def sheet   = file(params.input, checkIfExists: true)
    def sheetDir = sheet.parent
    def resolveRead = { p ->
        (p =~ /^([a-z]+:\/\/|\/)/) ? file(p, checkIfExists: true)
                                   : file("${sheetDir}/${p}", checkIfExists: true)
    }
    ch_reads = Channel
        .fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            def meta  = [ id: row.sample ]
            def reads = row.fastq_2 ? [ resolveRead(row.fastq_1), resolveRead(row.fastq_2) ]
                                    : [ resolveRead(row.fastq_1) ]
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
