process MULTIMM {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/sfglab/multimm:1.0':
        'docker.io/sfglab/multimm:1.0' }"
    
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
    tail -n +2 ${loops} > tmp_loops.bedpe
    
MultiMM \\
    --loops_path tmp_loops.bedpe \\
    --compartment_path ${compartment} \\
    --cpu_threads ${task.cpus} \\
    --out_path ./ \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MultiMM: 1.1.0.0.0
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
        MultiMM: 1.1.0.0.0
    END_VERSIONS
    """
}


