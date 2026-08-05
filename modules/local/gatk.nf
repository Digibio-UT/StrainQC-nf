/*
 * GATK variant-calling path. Ports backend/app/pipeline/wrappers.py:
 *   prepare_reference_fasta, run_bwa_mem, run_gatk_mark_duplicates,
 *   run_gatk_haplotype_caller, run_gatk_mutect2
 */

// gbk→fasta + samtools faidx + GATK dict + bwa-mem2 index
process PREPARE_REFERENCE {
    tag   { reference.baseName }
    label 'process_medium'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path reference

    output:
    tuple path("reference.fasta"), path("reference.fasta.*"), path("reference.dict"), emit: refset
    path "versions.yml", emit: versions

    script:
    """
    # Normalise reference to reference.fasta (GBK→FASTA if needed; skip self-copy when the
    # input is already named reference.fasta, which the staged symlink can be).
    case "${reference}" in
        *.gbk|*.gb|*.genbank) gbk_to_fasta.py "${reference}" reference.fasta ;;
        reference.fasta)      : ;;   # already the target name — use as-is
        *)                    cp "${reference}" reference.fasta ;;
    esac
    samtools faidx reference.fasta
    gatk --java-options "-Xmx${(int)(task.memory.toGiga()*0.8)}g" CreateSequenceDictionary -R reference.fasta -O reference.dict
    bwa-mem2 index reference.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK\\) v\\K[0-9.]+' || echo NA)
    END_VERSIONS
    """
    // NOTE: under -profile *,apptainer this needs samtools+bwa-mem2+gatk in one image.
    // Under hpc_modules the beforeScript loads all three. TODO: split for the container path.

    stub:
    """
    touch reference.fasta reference.fasta.fai reference.dict versions.yml
    """
}

process BWA_MEM2_ALIGN {
    tag   { meta.id }
    label 'process_high'
    publishDir "${params.outdir}/${meta.id}/gatk", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(meta), path(reads)
    tuple path(fasta), path(fai_and_index), path(dict)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1' \\
        ${fasta} ${reads} > aln.sam
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam aln.sam
    samtools index ${meta.id}.sorted.bam
    rm -f aln.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai versions.yml
    """
}

process GATK_MARKDUPLICATES {
    tag   { meta.id }
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), path("${meta.id}.dedup.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    gatk --java-options "-Xmx${(int)(task.memory.toGiga()*0.8)}g" MarkDuplicates \\
        -I ${bam} -O ${meta.id}.dedup.bam \\
        -M ${meta.id}.metrics.txt --REMOVE_DUPLICATES false --CREATE_INDEX true
    echo '"${task.process}": {gatk4: MarkDuplicates}' > versions.yml
    """

    stub:
    """
    touch ${meta.id}.dedup.bam ${meta.id}.dedup.bai versions.yml
    """
}

// HaplotypeCaller (consensus, needs ploidy) OR Mutect2 (polymorphism) — selected in the workflow.
process GATK_CALL {
    tag   { meta.id }
    label 'process_high'
    publishDir "${params.outdir}/${meta.id}/gatk", mode: 'copy', pattern: '*.vcf*'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai_and_index), path(dict)

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def minbq = "--min-base-quality-score ${params.gatk_min_base_quality}"
    def minmq = "--minimum-mapping-quality ${params.gatk_min_mapping_quality}"
    def maxr  = "--max-reads-per-alignment-start ${params.gatk_max_reads_per_align_start}"
    def soft  = params.gatk_dont_use_soft_clipped ? '--dont-use-soft-clipped-bases' : ''
    def xmx   = "-Xmx${(int)(task.memory.toGiga()*0.8)}g"
    if (params.polymorphism_mode)
        """
        gatk --java-options "${xmx}" Mutect2 -R ${fasta} -I ${bam} -O ${meta.id}.vcf \\
            ${minbq} ${minmq} ${maxr} ${soft} --callable-depth ${params.gatk_callable_depth}
        echo '"${task.process}": {gatk4: Mutect2}' > versions.yml
        """
    else
        """
        gatk --java-options "${xmx}" HaplotypeCaller -R ${fasta} -I ${bam} -O ${meta.id}.vcf \\
            --sample-ploidy ${params.ploidy} -ERC NONE \\
            ${minbq} ${minmq} ${maxr} ${soft} \\
            --standard-min-confidence-threshold-for-calling ${params.gatk_call_conf}
        echo '"${task.process}": {gatk4: HaplotypeCaller}' > versions.yml
        """

    stub:
    """
    touch ${meta.id}.vcf versions.yml
    """
}
