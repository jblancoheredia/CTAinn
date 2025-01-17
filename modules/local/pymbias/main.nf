process PYMBIAS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/pymbias:1.1.1' :
        'drgiovianco/pymbias:1.1.1' }"

    input:
    tuple val(meta), path(mbias)
    val(plot_type)
    val(plot_ax_x)
    val(plot_ax_y)

    output:
    tuple val(meta), path("*.tsv"),                 emit: cutoffs
    tuple val(meta), path("*.pdf"),  optional:true, emit: pdf
    tuple val(meta), path("*.svg"),  optional:true, emit: svg
    tuple val(meta), path("*.jpeg"), optional:true, emit: jpeg
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    PyMbias.py \\
        -p ${prefix} \\
        -s $plot_type \\
        -x $plot_ax_x \\
        -y $plot_ax_y \\
        $mbias

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymbias: \$(python PyMbias.py --version | sed 's/PyMbias v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "${plot_type}"
    """
    touch ${prefix}_mbias.tsv
    touch ${prefix}_mbias.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymbias: \$(python PyMbias.py --version | sed 's/PyMbias v//')
    END_VERSIONS
    """
}    
