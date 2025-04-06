
process COOLTOOLS_EIGSCIS {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cooltools:0.7.1--52a846d378a83a4f':
        'community.wave.seqera.io/library/cooltools:0.7.1--2b886db64cb0cd87' }"
    
    input:
    tuple val(meta), path(mcool)
    val(resolution)
    
    output:
    tuple val(meta), path("*_eigs_res_*.tsv"), emit: scores
    path "versions.yml"                 , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def version = "0.7.1"
    """
    cooltools eigs-cis \\
        ${mcool}::/resolutions/$resolution \\
        $args \\
        -o ${prefix}_eigs_res_${resolution}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools_eigs-cis: $version
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "0.7.1"
    
    """
    touch ${prefix}_eigs_res_${resolution}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools_eigs-cis: $version
    END_VERSIONS
    """
}
