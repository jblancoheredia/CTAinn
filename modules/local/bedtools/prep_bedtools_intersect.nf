process PREP_BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(intervals)
    val affix

    output:
    tuple val(meta), path("${meta.id}.${affix}.bam")    , emit: bam
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "${affix}.bam"
    def bamout = "${prefix}.${suffix}"
    """
    bedtools \\
        intersect \\
        -wa \\
        -u \\
        -abam $bam \\
        -b $intervals \\
        > ${bamout}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = "${affix}.bam"
    def bamout = "${prefix}.${suffix}"
    """
    touch ${bamout}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
