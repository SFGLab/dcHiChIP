
process JUICERTOOLS {
    tag "$meta.id"
    label 'process_medium'
    container "community.wave.seqera.io/library/openjdk:23.0.2--2fd1f5d679ee38ac"
    
    input:
    tuple val(meta), path(hic_input)
    
    output:
    tuple val(meta),       path("*.hic"), emit: hic
    path "versions.yml"                       , emit: versions
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    juicer_tools pre \\
    ${hic_input} \\
    ${prefix}.hic \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_tools: 2.20.00
    END_VERSIONS
    """
}
