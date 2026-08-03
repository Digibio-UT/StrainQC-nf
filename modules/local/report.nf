/*
 * Final report process — the single point that emits results.json.
 * Runs bin/emit_results.py, the port of backend/app/pipeline/parsers.py + manager.py
 * assembly. See docs/API_CONTRACT.md for the exact schema. (`manta_vcf` intentionally absent.)
 *
 * Optional inputs use ONE empty placeholder file (assets/NO_FILE): emit_results.py treats
 * any empty / NO_*-named path as "that step didn't run" (matches manager.py's os.path.exists
 * guards), so skipped steps (no Kraken DB, breseq path, etc.) just pass the placeholder.
 */
process EMIT_RESULTS {
    tag   { meta.id }
    label 'process_low'
    container params.containers.report
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta),
          path(variant_src),                                // breseq 'output/' dir OR gatk .vcf
          path(reference),                                  // reference (name → results.reference_file)
          path(coverage,     stageAs: 'coverage.tsv'),
          path(flagstat,     stageAs: 'flagstat.txt'),
          path(kraken,       stageAs: 'kraken.report'),
          path(mash,         stageAs: 'mash_screen.tab'),
          path(mash_unmap,   stageAs: 'mash_unmapped.tab'),
          path(blast,        stageAs: 'blast_results.tsv'),
          path(evidence,     stageAs: 'evidence'),          // dir of per-variant JSON
          path(unmapped)                                    // 0..n FASTQ (names → unmapped_files)

    output:
    tuple val(meta), path('results.json'), emit: results
    path 'high_confidence_snps.csv', optional: true
    path 'coverage.tsv',             optional: true         // re-published for the /coverage endpoint

    script:
    def variant_arg  = params.variant_caller == 'gatk' ? "--vcf ${variant_src}"
                                                        : "--breseq-dir ${variant_src}"
    def unmapped_arg = unmapped ? "--unmapped ${unmapped}" : ''
    """
    emit_results.py \\
        --job-id ${meta.id} \\
        --caller ${params.variant_caller} \\
        --polymorphism ${params.polymorphism_mode} \\
        --reference ${reference} --reference-gbk ${reference} \\
        ${variant_arg} \\
        --coverage coverage.tsv \\
        --flagstat flagstat.txt \\
        --kraken kraken.report \\
        --mash mash_screen.tab --mash-unmapped mash_unmapped.tab \\
        --blast blast_results.tsv \\
        --evidence-dir evidence \\
        ${unmapped_arg} \\
        --out results.json
    """

    stub:
    """
    echo '{"job_id":"${meta.id}","status":"completed","variant_caller":"${params.variant_caller}","mutations":[]}' > results.json
    """
}
