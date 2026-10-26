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
    if (!(params.variant_caller in ['gatk', 'breseq']))
        error "--variant_caller must be 'gatk' or 'breseq' (got '${params.variant_caller}')"

    // Read the samplesheet up front so it can be validated before any task runs.
    def sheet    = file(params.input, checkIfExists: true)
    def sheetDir = sheet.parent
    def rows     = sheet.splitCsv(header: true)
    if (!rows || !rows[0].containsKey('sample') || !rows[0].containsKey('fastq_1'))
        error "Samplesheet ${sheet.name} needs a header with at least: sample,fastq_1[,fastq_2]"
    def ids   = rows.collect { it.sample }
    def dupes = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (dupes)
        error "Duplicate sample ids in ${sheet.name}: ${dupes.join(', ')} — ids key every output folder and join"

    // CNV runs only when the named parent is actually in the samplesheet; a typo would
    // otherwise skip CNV silently while the run still succeeds.
    if (params.parent_sample) {
        if (!(params.parent_sample in ids))
            error "--parent_sample '${params.parent_sample}' is not a sample in ${sheet.name} (samples: ${ids.join(', ')})"
        if (params.variant_caller != 'gatk')
            log.warn "--parent_sample is set but CNV calling runs only on the GATK path; it is ignored for '${params.variant_caller}'"
    }

    // Parse rows → [ meta, [reads] ]. Read paths may be absolute, a URL, or relative to the
    // samplesheet's own directory (keeps bundled test fixtures portable).
    def resolveRead = { p ->
        (p =~ /^([a-z]+:\/\/|\/)/) ? file(p, checkIfExists: true)
                                   : file("${sheetDir}/${p}", checkIfExists: true)
    }
    ch_reads = Channel
        .fromList(rows)
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
