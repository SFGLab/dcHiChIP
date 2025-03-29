
process MULTIMM {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_multimm:8c581e45aa852599':
        'community.wave.seqera.io/library/pip_multimm:a503e99c16691d6a' }"
    
    input:
    tuple val(meta), path(fastq1), path(fastq2)
    tuple val(meta2), path(peaks)
    tuple val(meta3), path(bwa_index)
    path (features)
    path (reqiured_scripts)

    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    tuple val(meta), path("*.hic.input"), emit: hic
    path "versions.yml"                 , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    MultiMM -c config.ini ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
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
        MAPS: 1.1.0
    END_VERSIONS
    """
}
