/*
 * samtools QC + evidence pre-compute. Ports wrappers.py:
 *   run_samtools_depth, run_samtools_flagstat, run_samtools_extract_unmapped,
 *   extract_variant_read_evidence
 *
 * KEY re-platform decision (open strategy Q2): coverage downsampling and per-variant
 * read-evidence JSON are PRE-COMPUTED here in-pipeline, so the API only serves files
 * instead of re-reading BAMs on demand (expensive on cloud object storage).
 */

process SAMTOOLS_DEPTH {
    tag   { meta.id }
    label 'process_low'
    container params.containers.samtools

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("coverage.tsv"), emit: coverage

    script:
    """
    samtools depth -a ${bam} > coverage.tsv
    """

    stub:
    "touch coverage.tsv"
}

process SAMTOOLS_FLAGSTAT {
    tag   { meta.id }
    label 'process_low'
    container params.containers.samtools

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("flagstat.txt"), emit: flagstat

    script:
    "samtools flagstat ${bam} > flagstat.txt"

    stub:
    "touch flagstat.txt"
}

process SAMTOOLS_UNMAPPED {
    tag   { meta.id }
    label 'process_low'
    container params.containers.samtools
    publishDir "${params.outdir}/${meta.id}/gatk", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("unmapped.fastq"), emit: unmapped, optional: true

    script:
    """
    samtools view -b -f 4 -o unmapped.bam ${bam}
    samtools fastq unmapped.bam -o unmapped.fastq || true
    [ -s unmapped.fastq ] || rm -f unmapped.fastq
    """

    stub:
    "touch unmapped.fastq"
}

// Pre-compute per-variant read-evidence JSON (the /evidence endpoint just serves these).
process VARIANT_EVIDENCE {
    tag   { meta.id }
    label 'process_medium'
    container params.containers.report
    publishDir "${params.outdir}/${meta.id}/gatk", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    tuple path(fasta), path(fai_and_index), path(dict)

    output:
    tuple val(meta), path("evidence"), emit: evidence

    script:
    """
    mkdir -p evidence
    # TODO: port extract_variant_read_evidence() from wrappers.py into bin/extract_evidence.py
    #   inputs:  ${bam}, ${fasta}, mutations (from ${vcf}), window=${params.evidence_window}, max_reads=${params.evidence_max_reads}
    #   output:  evidence/<chrom>_<start>_<end>_<n>.json  (schema unchanged from current app)
    extract_evidence.py --bam ${bam} --ref ${fasta} --vcf ${vcf} \\
        --window ${params.evidence_window} --max-reads ${params.evidence_max_reads} \\
        --out evidence || true
    """

    stub:
    "mkdir -p evidence"
}
