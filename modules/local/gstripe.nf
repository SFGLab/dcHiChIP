process GSTRIPE {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4eaae53e6b938c3212eb34fc89f4c1d514511a3e05bb3f91f5b60a41e6870f9b/data' :
        'community.wave.seqera.io/library/pip_gstripe:6d10173fe52641bf'}"

    input:
    tuple val(meta) , path(bedpe)

    output:
    tuple val(meta), path("*gstripes_raw.tsv"),  emit: gstripe
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 -m gstripe.gstripe \\
        ${bedpe} \\
        ./ \\
        --max_workers=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${bedpe}.gstripes_raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
