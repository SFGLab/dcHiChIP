process GSTRIPE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioframe=0.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        '' }"

    input:
    tuple val(meta) , path(bedpe)
    
    output:
    tuple val(meta), path("*_gstripe.txt"),  emit: gstripe
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    python -m gstripe.gstripe \\
        ${bedpe} \\
        ${prefix}_gstripe.txt \\
        --max_workers=${task.cpus}

    with open ("versions.yml", "w") as version_file:
	   version_file.write("\\"${task.process}\\":\\n    python: {}\\n".format(sys.version.split()[0].strip()))

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}_gstripe.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
