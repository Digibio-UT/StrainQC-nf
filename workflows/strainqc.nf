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
include { SAMTOOLS_DEPTH; SAMTOOLS_FLAGSTAT; SAMTOOLS_UNMAPPED;
          VARIANT_EVIDENCE; BRESEQ_COVERAGE }            from '../modules/local/qc.nf'
include { KRAKEN2;
          MASH_SCREEN as MASH_RAW;
          MASH_SCREEN as MASH_UNMAPPED;
          BLASTN_REMOTE }                               from '../modules/local/contamination.nf'
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

        // Contamination (matches old app): Kraken2 on UNMAPPED → contaminants;
        // MASH on RAW reads → mash_results (identifies dominant + contaminant strains);
        // MASH on UNMAPPED → mash_unmapped_results. Each skipped if its DB param is unset.
        def kdb = params.kraken2_db ? file(params.kraken2_db) : []
        def mdb = params.mash_db    ? file(params.mash_db)    : []
        KRAKEN2(SAMTOOLS_UNMAPPED.out.unmapped, kdb)
        MASH_RAW(ch_reads, mdb)
        MASH_UNMAPPED(SAMTOOLS_UNMAPPED.out.unmapped, mdb)

        // Assemble the 11-slot EMIT_RESULTS tuple keyed by meta, filling gaps with NF.
        // (reference is one value for all samples → appended in the map, not joined.)
        ch_report = GATK_CALL.out.vcf                                       // [meta, vcf]
            .join(SAMTOOLS_DEPTH.out.coverage,       remainder: true)
            .join(SAMTOOLS_FLAGSTAT.out.flagstat,    remainder: true)
            .join(KRAKEN2.out.report,                remainder: true)
            .join(MASH_RAW.out.screen,               remainder: true)
            .join(MASH_UNMAPPED.out.screen,          remainder: true)
            .join(VARIANT_EVIDENCE.out.evidence,     remainder: true)
            .join(SAMTOOLS_UNMAPPED.out.unmapped,    remainder: true)
            .map { meta, vcf, cov, flag, krak, mashraw, mashunmap, evid, unmap ->
                tuple(meta, vcf, ch_reference,
                      orNF(cov), orNF(flag), orNF(krak),
                      orNF(mashraw), orNF(mashunmap),
                      NF /*blast fallback: TODO wire BLASTN_REMOTE when no kraken DB*/,
                      orNF(evid), orNF(unmap))
            }

    } else {   // ── breseq path ───────────────────────────────────────
        BRESEQ(ch_reads, ch_reference)

        // Coverage from breseq's own BAM (so results.json gets coverage_file → the frontend
        // renders the coverage plot + genome browser). read_stats/evidence come from the
        // breseq output/ tree, which emit_results.py reads directly from breseq-dir.
        BRESEQ_COVERAGE(BRESEQ.out.outdir)

        // contamination on RAW reads (breseq path; breseq's own unmapped extraction TODO)
        KRAKEN2(ch_reads, params.kraken2_db ? file(params.kraken2_db) : [])
        MASH_RAW(ch_reads, params.mash_db ? file(params.mash_db) : [])

        ch_report = BRESEQ.out.outdir                                       // [meta, outdir]
            .join(BRESEQ_COVERAGE.out.coverage, remainder: true)
            .join(KRAKEN2.out.report,           remainder: true)
            .join(MASH_RAW.out.screen,          remainder: true)
            .map { meta, outdir, cov, krak, mash ->
                tuple(meta, outdir, ch_reference,
                      orNF(cov), NF /*flagstat: breseq read_stats come from summary.html*/,
                      orNF(krak), orNF(mash),
                      NF /*mash_unmapped*/, NF /*blast*/, NF /*evidence*/, NF /*unmapped*/)
            }
    }

    EMIT_RESULTS(ch_report)

    emit:
    results = EMIT_RESULTS.out.results
}
