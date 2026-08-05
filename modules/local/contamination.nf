/*
 * Contamination screening. Ports wrappers.py:
 *   run_kraken2, run_mash_screen, run_blast_sample / run_blastn
 * Any of these may be skipped depending on which reference DBs are staged at the site.
 */

process KRAKEN2 {
    tag   { meta.id }
    label 'process_high'
    publishDir "${params.outdir}/${meta.id}/contamination", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path("kraken.report"), emit: report

    when:
    db

    script:
    """
    kraken2 --db ${db} --output kraken.output --report kraken.report \\
        --threads ${task.cpus} ${reads}
    """

    stub:
    "touch kraken.report"
}

process MASH_SCREEN {
    tag   { meta.id }
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}/contamination", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path  msh

    output:
    tuple val(meta), path("mash_screen.tab"), emit: screen

    when:
    msh

    script:
    """
    mash screen -p ${task.cpus} -w ${msh} ${reads} > mash_screen.tab
    """

    stub:
    "touch mash_screen.tab"
}

process BLASTN_REMOTE {
    tag   { meta.id }
    label 'process_network'          // network-bound; 1 slot, long wall clock
    publishDir "${params.outdir}/${meta.id}/contamination", mode: 'copy'

    input:
    tuple val(meta), path(query_fasta)

    output:
    tuple val(meta), path("blast_results.tsv"), emit: hits

    when:
    params.run_blast

    script:
    """
    blastn -db nt -query ${query_fasta} -out blast_results.tsv -remote \\
        -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' \\
        -max_target_seqs 1 -max_hsps 1
    """

    stub:
    "touch blast_results.tsv"
}
