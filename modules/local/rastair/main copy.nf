process RASTAIR {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/rastair:0.8.0' :
        'drgiovianco/rastair:0.8.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(metdir)

    output:
    tuple val(meta), path("*.mods")         , emit: mods
    tuple val(meta), path("*.mods.summary") , emit: summary
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = "${metdir}/hg38_pUC19_lambda_unmod2kb.fa"

    """
    rastair \\
    call \\
    --nOT 0,0,0,0 \\
    --nOB 0,0,0,0 \\
    --fasta-file ${reference} \\
    $bam \\
    >> ${prefix}_rastair_output.mods

    rastair_summarize.py \\
        --rastair_in  ${prefix}_rastair_output.mods \\
        --summary_out ${prefix}_rastair_output.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """
}
