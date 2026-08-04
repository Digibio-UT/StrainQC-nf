/*
 * breseq — consensus & polymorphism variant calling for haploid microbial genomes.
 * Ports backend/app/pipeline/wrappers.py :: run_breseq
 * Emits the breseq output tree (HTML + output.gd + data/) consumed by the report step.
 */
process BRESEQ {
    tag   { meta.id }
    label 'process_high'
    container params.containers.breseq
    // Publish the breseq tree at results/<id>/breseq/ so its static evidence HTML resolves
    // at /results_files/<id>/breseq/output/evidence/… (the URL the frontend builds). The
    // breseq -o dir is therefore named `breseq` and published into the job root.
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  reference

    output:
    tuple val(meta), path("breseq"), emit: outdir
    tuple val(meta), path("breseq/output/output.gd"), emit: gd, optional: true
    path  "versions.yml", emit: versions

    script:
    def poly   = params.polymorphism_mode ? '-p' : ''
    def limit  = params.breseq_limit_fold_coverage  > 0 ? "--limit-fold-coverage ${params.breseq_limit_fold_coverage}" : ''
    def baseq  = params.breseq_base_quality_cutoff != 3 ? "--base-quality-cutoff ${params.breseq_base_quality_cutoff}" : ''
    def minmq  = params.breseq_min_mapping_quality != 0 ? "--minimum-mapping-quality ${params.breseq_min_mapping_quality}" : ''
    def polyfq = (params.polymorphism_mode && params.breseq_polymorphism_frequency_cutoff != 0.05) ? "--polymorphism-frequency-cutoff ${params.breseq_polymorphism_frequency_cutoff}" : ''
    def nojunc = params.breseq_no_junction_prediction ? '--no-junction-prediction' : ''
    def covstr = params.breseq_consensus_min_cov_each_strand > 0 ? "--consensus-minimum-total-coverage-each-strand ${params.breseq_consensus_min_cov_each_strand}" : ''
    def target = params.breseq_targeted_sequencing ? '--targeted-sequencing' : ''
    """
    mkdir -p breseq/output/evidence
    breseq -r ${reference} ${reads} -o breseq -j ${task.cpus} \\
        ${poly} ${limit} ${baseq} ${minmq} ${polyfq} ${nojunc} ${covstr} ${target}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        breseq: \$(breseq --version 2>&1 | head -n1 | sed 's/^breseq //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p breseq/output/evidence
    touch breseq/output/output.gd versions.yml
    """
}
