process CCD_CALLER {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/523d0e75031a5d95b56a4daca9d70fdd3e857c38c7e403d1f1ab04141d3fa539/data' :
        'community.wave.seqera.io/library/pip_numpy_pandas:9b0b2cb6e7b996bb' }"

    input:
    tuple val(meta), path(bedpe)

    output:
    tuple val(meta), path("*.ccds.bed"),  emit: ccds
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ccd_calling.py \\
    ${bedpe} \\
    ${args} \\
    -o ${prefix}.ccds.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.ccds.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
