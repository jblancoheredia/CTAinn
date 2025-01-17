process BWA_METH {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/bwameth:0.2.7' :
        'blancojmskcc/bwameth:0.2.7' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bwameth.py \\
        $args \\
        -t $task.cpus \\
        --reference $index/$fasta \\
        $reads \\
        | samtools sort $args2 -@ ${task.cpus-1} -o ${prefix}_sorted.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwameth: \$(bwameth.py --version | cut -f2 -d" ")
    END_VERSIONS
    """
}
