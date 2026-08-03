/*
 * Final report process — the single point that emits results.json.
 * Ports backend/app/pipeline/parsers.py (breseq HTML/GD, kraken report, mash, VCF).
 *
 * CONTRACT: results.json MUST match docs/API_CONTRACT.md exactly on every backend,
 * so the existing React frontend + FastAPI serving layer work unchanged.
 * (Note: `manta_vcf` was removed from the contract in the working app — do NOT emit it.)
 */
process EMIT_RESULTS {
    tag   { meta.id }
    label 'process_low'
    container params.containers.report
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta),
          path(variant_source),      // breseq output/ dir  OR  gatk vcf
          path(coverage,   stageAs: 'coverage.tsv'),
          path(flagstat,   stageAs: 'flagstat.txt'),
          path(kraken,     stageAs: 'kraken.report'),
          path(mash,       stageAs: 'mash_screen.tab'),
          path(evidence)               // evidence/ dir (may be empty)

    output:
    tuple val(meta), path("results.json"), emit: results

    script:
    """
    emit_results.py \\
        --job-id       ${meta.id} \\
        --caller       ${params.variant_caller} \\
        --polymorphism ${params.polymorphism_mode} \\
        --variant-src  ${variant_source} \\
        --coverage     coverage.tsv \\
        --flagstat     flagstat.txt \\
        --kraken       kraken.report \\
        --mash         mash_screen.tab \\
        --evidence-dir ${evidence} \\
        --coverage-resolution ${params.coverage_resolution} \\
        --out results.json
    """

    stub:
    """
    echo '{"job_id":"${meta.id}","status":"completed","variant_caller":"${params.variant_caller}","mutations":[]}' > results.json
    """
}
