process BWAMEM2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' :
        'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' }"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta1), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    val(sort_bam)

    output:
    tuple val(meta), path("*.bam") , path("*.bai"), emit: bam
    tuple val(meta), path("*.csi") ,                emit: csi , optional:true
    tuple val(meta), path("*.sam") ,                emit: sam , optional:true
    tuple val(meta), path("*.cram"),                emit: cram, optional:true
    tuple val(meta), path("*.crai"),                emit: crai, optional:true
    path  "versions.yml"           ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    cat ${conss_r1} ${split_r1} > ${prefix}_R1.fq.gz
    cat ${conss_r2} ${split_r2} > ${prefix}_R2.fq.gz

    bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        ${reads} \\
        | samtools $samtools_command $args2 -@ $task.cpus ${reference} -o ${prefix}.${extension} -

    samtools index -@ $task.cpus -o ${prefix}.${extension}.bai ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:

    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.${extension}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
