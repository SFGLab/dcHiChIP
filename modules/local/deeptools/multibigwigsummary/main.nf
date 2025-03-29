process DEEPTOOLS_MUTIBIGWIGSUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(big_wigs)
    val(labels)

    output:
    tuple val(meta), path("*.npz"), emit: npz
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    multiBigwigSummary bins \\
        $args \\
        --bwfiles ${big_wigs} \\
        --outFileName ${prefix}_results.npz \\
        --numberOfProcessors ${task.cpus} \\
        --labels ${labels.join(" ")}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed -e "s/multiBigwigSummary //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.results.npz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(multiBigwigSummary --version | sed -e "s/multiBigwigSummary //g")
    END_VERSIONS
    """
}
