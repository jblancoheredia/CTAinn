process RASTAIR {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/rastair:0.8.0' :
        'drgiovianco/rastair:0.8.0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    val(cutoffs)

    output:
    tuple val(meta), path("*.mods")         , emit: mods
    tuple val(meta), path("*.mods.summary") , emit: summary
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rastair \\
        call \\
        --nOT $cutoffs.not \\
        --nOB $cutoffs.nob \\
        --fasta-file ${fasta} \\
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
    touch ${prefix}_rastair_output.mods
    touch ${prefix}_rastair_output.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """
}

process RASTAIR_SUMMARY {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/rastair:0.8.0' :
        'drgiovianco/rastair:0.8.0' }"

    input:
    tuple val(meta), path(mods)

    output:
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("*.targeted.mods.summary"), emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rastair_summarize.py \\
        --rastair_in  $mods \\
        --summary_out ${prefix}_rastair_output.targeted.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rastair_output.targeted.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """
}

process RASTAIR_PERREAD {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/rastair:0.8.0' :
        'drgiovianco/rastair:0.8.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)

    output:
    path "versions.yml"                             , emit: versions
    tuple val(meta), path("*_perread.mods.summary"), emit: summary
    tuple val(meta), path("*_perread.mods")         , emit: mods

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rastair \\
        per-read \\
        --fasta-file ${fasta} \\
        $bam \\
        >> ${prefix}_rastair_output_perread.mods

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rastair_output_perread.mods
    touch ${prefix}_rastair_output_perread.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """
}

process RASTAIR_MBIAS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/rastair:0.8.0' :
        'drgiovianco/rastair:0.8.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)

    output:
    path "versions.yml"             , emit: versions
    tuple val(meta), path("*.mbias"), emit: mbias

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rastair \\
        mbias \\
        --fasta-file ${fasta} \\
        $bam \\
        >> ${prefix}_rastair_output.mbias



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rastair_output.mbias

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version |& sed '1!d ; s/rastair //')
    END_VERSIONS
    """
}