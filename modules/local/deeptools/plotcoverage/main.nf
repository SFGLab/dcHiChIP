process DEEPTOOLS_PLOTCOVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bams), path (index)
    val(labels)
    output:
    tuple val(meta), path("*.png"), emit: coverage
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotCoverage \\
    --bamfiles ${bams} \\
    --plotFile ${prefix}_coverage.png \\
    --numberOfProcessors ${task.cpus} \\
    --labels ${labels.join(" ")} \\
    $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotCoverage --version | sed -e "s/plotCoverage //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_coverage.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotCoverage --version | sed -e "s/plotCoverage //g")
    END_VERSIONS
    """
}
