process MULTIMM {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_multimm:8c581e45aa852599':
        'community.wave.seqera.io/library/pip_multimm:a503e99c16691d6a' }"
    
    input:
    tuple val(meta), path(loops)
    tuple val(meta2), path(compartment)
    
    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    tuple val(meta), path("*.hic.input"), emit: hic
    path "versions.yml"                 , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    MultiMM \\
    --loops_path ${loops} \\
    --compartment_path ${compartment} \\
    --cpu_threads ${task.cpus} \\
    --out_path ./ \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MultiMM: 1.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.bedpe
    touch ${prefix}.hic.input
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MultiMM: 1.1.0
    END_VERSIONS
    """
}
