/*
 * Main StrainQC workflow. Two variant-calling paths (breseq | gatk) converge on a
 * single EMIT_RESULTS process that produces the frozen results.json contract.
 */

include { BRESEQ }                                       from '../modules/local/breseq.nf'
include { PREPARE_REFERENCE; BWA_MEM2_ALIGN;
          GATK_MARKDUPLICATES; GATK_CALL }               from '../modules/local/gatk.nf'
include { SAMTOOLS_DEPTH; SAMTOOLS_FLAGSTAT;
          SAMTOOLS_UNMAPPED; VARIANT_EVIDENCE }          from '../modules/local/qc.nf'
include { KRAKEN2; MASH_SCREEN; BLASTN_REMOTE }          from '../modules/local/contamination.nf'
include { EMIT_RESULTS }                                 from '../modules/local/report.nf'

workflow STRAINQC {

    take:
    ch_reads          // [ meta, [reads] ]
    ch_reference      // path

    main:
    def PLACEHOLDER = file("${projectDir}/assets/NO_FILE")   // for optional report inputs

    if (params.variant_caller == 'gatk') {
        PREPARE_REFERENCE(ch_reference)
        def refset = PREPARE_REFERENCE.out.refset

        BWA_MEM2_ALIGN(ch_reads, refset)
        GATK_MARKDUPLICATES(BWA_MEM2_ALIGN.out.bam)
        def bam = GATK_MARKDUPLICATES.out.bam
        GATK_CALL(bam, refset)

        SAMTOOLS_DEPTH(bam)
        SAMTOOLS_FLAGSTAT(bam)
        SAMTOOLS_UNMAPPED(bam)

        // evidence needs bam + vcf together
        def bam_vcf = bam.join(GATK_CALL.out.vcf)          // [meta, bam, bai, vcf]
        VARIANT_EVIDENCE(bam_vcf, refset)

        // contamination runs on UNMAPPED reads in the GATK path
        KRAKEN2(SAMTOOLS_UNMAPPED.out.unmapped.map { m, f -> [m, f] }, params.kraken2_db ? file(params.kraken2_db) : [])
        MASH_SCREEN(SAMTOOLS_UNMAPPED.out.unmapped, params.mash_db ? file(params.mash_db) : [])

        // assemble report inputs (variant_source = vcf)
        ch_report = GATK_CALL.out.vcf
            .join(SAMTOOLS_DEPTH.out.coverage)
            .join(SAMTOOLS_FLAGSTAT.out.flagstat)
            .join(KRAKEN2.out.report,   remainder: true)
            .join(MASH_SCREEN.out.screen, remainder: true)
            .join(VARIANT_EVIDENCE.out.evidence, remainder: true)
            .map { it.collect { v -> v ?: PLACEHOLDER } }

    } else {   // ── breseq path ───────────────────────────────────────
        BRESEQ(ch_reads, ch_reference)

        // contamination runs on the RAW reads in the breseq path
        KRAKEN2(ch_reads, params.kraken2_db ? file(params.kraken2_db) : [])
        MASH_SCREEN(ch_reads, params.mash_db ? file(params.mash_db) : [])

        // breseq produces its own coverage/flagstat/evidence inside output/ — pass placeholders
        // for the samtools-derived files; the report parser reads them from the breseq tree.
        ch_report = BRESEQ.out.outdir
            .map { m, d -> [m, d, PLACEHOLDER, PLACEHOLDER] }
            .join(KRAKEN2.out.report,     remainder: true)
            .join(MASH_SCREEN.out.screen, remainder: true)
            .map { it + [PLACEHOLDER] }   // evidence dir slot
            .map { it.collect { v -> v ?: PLACEHOLDER } }
    }

    EMIT_RESULTS(ch_report)

    emit:
    results = EMIT_RESULTS.out.results
}
