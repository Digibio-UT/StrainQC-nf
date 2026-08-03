/*
 * Main StrainQC workflow. Two variant-calling paths (breseq | gatk) converge on a
 * single EMIT_RESULTS process that produces the frozen results.json contract.
 *
 * ⚠️  The EMIT_RESULTS input is an 11-slot tuple (see modules/local/report.nf). Absent
 *     optional slots are filled with the empty placeholder assets/NO_FILE, which
 *     emit_results.py reads as "step skipped". This channel wiring has NOT yet been
 *     validated with a real Nextflow run — do `nextflow run . -profile test,docker -stub`
 *     as the next step (see README Roadmap).
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
    // Empty placeholder for optional EMIT_RESULTS slots (emit_results.py reads empty → skipped).
    def NF = file("${projectDir}/assets/NO_FILE")

    // Helper: coalesce a joined channel value (may be null from remainder joins) to NF.
    def orNF = { v -> v ?: NF }

    if (params.variant_caller == 'gatk') {
        PREPARE_REFERENCE(ch_reference)
        def refset = PREPARE_REFERENCE.out.refset

        BWA_MEM2_ALIGN(ch_reads, refset)
        GATK_MARKDUPLICATES(BWA_MEM2_ALIGN.out.bam)
        def bam = GATK_MARKDUPLICATES.out.bam                 // [meta, bam, bai]
        GATK_CALL(bam, refset)

        SAMTOOLS_DEPTH(bam)
        SAMTOOLS_FLAGSTAT(bam)
        SAMTOOLS_UNMAPPED(bam)

        def bam_vcf = bam.join(GATK_CALL.out.vcf)             // [meta, bam, bai, vcf]
        VARIANT_EVIDENCE(bam_vcf, refset)

        // contamination on UNMAPPED reads (GATK path)
        KRAKEN2(SAMTOOLS_UNMAPPED.out.unmapped, params.kraken2_db ? file(params.kraken2_db) : [])
        MASH_SCREEN(SAMTOOLS_UNMAPPED.out.unmapped, params.mash_db ? file(params.mash_db) : [])

        // Assemble the 11-slot EMIT_RESULTS tuple keyed by meta, filling gaps with NF.
        ch_report = GATK_CALL.out.vcf                                       // [meta, vcf]
            .combine(ch_reference)                                          // + reference
            .join(SAMTOOLS_DEPTH.out.coverage,       remainder: true)
            .join(SAMTOOLS_FLAGSTAT.out.flagstat,    remainder: true)
            .join(KRAKEN2.out.report,                remainder: true)
            .join(MASH_SCREEN.out.screen,            remainder: true)
            .join(VARIANT_EVIDENCE.out.evidence,     remainder: true)
            .join(SAMTOOLS_UNMAPPED.out.unmapped,    remainder: true)
            .map { meta, vcf, ref, cov, flag, krak, mash, evid, unmap ->
                tuple(meta, vcf, ref,
                      orNF(cov), orNF(flag), orNF(krak), orNF(mash),
                      NF /*mash_unmapped: TODO wire second MASH on unmapped*/,
                      NF /*blast fallback: TODO wire BLASTN_REMOTE when no kraken DB*/,
                      orNF(evid), orNF(unmap))
            }

    } else {   // ── breseq path ───────────────────────────────────────
        BRESEQ(ch_reads, ch_reference)

        // contamination on RAW reads (breseq path)
        KRAKEN2(ch_reads, params.kraken2_db ? file(params.kraken2_db) : [])
        MASH_SCREEN(ch_reads, params.mash_db ? file(params.mash_db) : [])

        // breseq produces coverage/read-stats/evidence inside its own output/ tree, so the
        // samtools-derived slots are placeholders; emit_results.py reads them from breseq-dir.
        ch_report = BRESEQ.out.outdir                                       // [meta, outdir]
            .combine(ch_reference)                                          // + reference
            .join(KRAKEN2.out.report,     remainder: true)
            .join(MASH_SCREEN.out.screen, remainder: true)
            .map { meta, outdir, ref, krak, mash ->
                tuple(meta, outdir, ref,
                      NF /*coverage*/, NF /*flagstat*/,
                      orNF(krak), orNF(mash),
                      NF /*mash_unmapped*/, NF /*blast*/, NF /*evidence*/, NF /*unmapped*/)
            }
    }

    EMIT_RESULTS(ch_report)

    emit:
    results = EMIT_RESULTS.out.results
}
