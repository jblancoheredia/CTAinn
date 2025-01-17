process ASTAIR {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://drgiovianco/astair:3.3.2' :
        'drgiovianco/astair:3.3.2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(ref_dir)
    val(read_length)

    output:
    tuple val(meta), path("*.mods")  , emit: mods
    tuple val(meta), path("*.stats") , emit: stats

    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_astair_output
    rm ${fasta} ${fasta_fai}
    cp ${ref_dir}/${fasta} .
    cp ${ref_dir}/${fasta_fai} .

    astair \\
        call \\
        --N_threads $task.cpus \\
        --no_information 0 \\
        --skip_clip_overlap true \\
        --start_clip 7 \\
        --end_clip 7 \\
        -i $bam \\
        -f ${fasta} \\
        --context CpG \\
        -d ./${prefix}_astair_output

    astair \\
        mbias \\
        -i $bam \\
        -f ${fasta} \\
        --read_length $read_length \\
        --plot \\
        -d ./${prefix}_astair_output

    mv ${prefix}_astair_output/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        astair: \$(astair --help | grep "Version" | sed 's/^[ \t]*Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_astair_output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        astair: \$(astair --help | grep "Version" | sed 's/^[ \t]*Version: //')
    END_VERSIONS
    """
}    
