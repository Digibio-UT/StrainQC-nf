/*
 * Copy-number / aneuploidy calling with CNVkit (WGS mode).
 *
 * Cohort step, GATK path only. Mirrors the Ottilie et al. (2022) approach of normalising
 * evolved-clone coverage against the un-evolved PARENT: the parent BAM builds a CNVkit
 * reference (CNVKIT_REFERENCE), then every non-parent sample is called against it
 * (CNVKIT_CALL). bin/summarize_cnv.py rolls the per-segment log2 ratios up to a per-chromosome
 * mean copy number and flags whole-chromosome aneuploidy → cnv.json (an optional EMIT_RESULTS
 * slot; absent for the parent and for any run without --parent_sample).
 *
 * Enabled only when params.parent_sample is set, so single-sample UI runs and the synthetic
 * `test` profile are unaffected.
 */

// Build a CNVkit reference (.cnn) from the parent/ancestor BAM (the "normal").
process CNVKIT_REFERENCE {
    tag   { meta.id }
    label 'process_medium'
    publishDir "${params.outdir}/cnv", mode: 'copy', pattern: 'reference.cnn'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai_and_index), path(dict)

    output:
    path "reference.cnn", emit: reference

    script:
    """
    samtools index ${bam}
    cnvkit.py batch \\
        --normal ${bam} \\
        --method wgs \\
        --fasta ${fasta} \\
        --output-reference reference.cnn \\
        --processes ${task.cpus} \\
        --output-dir .
    """

    stub:
    "touch reference.cnn"
}

// Call CNVs for one evolved sample against the parent-derived reference, then summarise
// to per-chromosome copy number + aneuploidy calls.
process CNVKIT_CALL {
    tag   { meta.id }
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}/cnv", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference
    tuple path(fasta), path(fai_and_index), path(dict)

    output:
    tuple val(meta), path("cnv.json"),      emit: cnv
    tuple val(meta), path("${meta.id}.cns"), emit: segments, optional: true

    script:
    // cnvkit names outputs after the BAM basename: <base>.cns is the segmented call file we
    // want; <base>.bintest.cns is a separate per-bin test (must NOT be picked up instead).
    // summarize_cnv.py calls aneuploidy by per-chromosome CN relative to the genome median
    // (robust to mosaic gains, e.g. Ottilie CBR110 chr I ~1.3x).
    def base = bam.baseName
    """
    samtools index ${bam}
    cnvkit.py batch ${bam} \\
        --reference ${reference} \\
        --method wgs \\
        --processes ${task.cpus} \\
        --output-dir .

    cp ${base}.cns ${meta.id}.cns

    # Per-chromosome copy number + whole-chromosome aneuploidy flags → cnv.json.
    summarize_cnv.py \\
        --cns ${meta.id}.cns \\
        --fai ${fasta}.fai \\
        --sample ${meta.id} \\
        --baseline-ploidy ${params.ploidy} \\
        --out cnv.json
    """

    stub:
    """
    touch ${meta.id}.cns
    echo '{"cnv_tool":"cnvkit","baseline_ploidy":${params.ploidy},"chromosome_copy_number":[],"aneuploidy":[]}' > cnv.json
    """
}
