process COOLTOOLS_BED_INVERT {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/pip_pandas_scipy:4f0dbe7431dea393' :
        'community.wave.seqera.io/library/pip_pandas_scipy:4f0dbe7431dea393' }"

    input:
    tuple val(meta) , path(bed_in)

    output:
    tuple val(meta), path("*_compartments.bed"),  emit: compartments
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in sfglab/dchichip/bin/
    """
    cooltools_bed_invert.py \\
        --path_comp ${bed_in} \\
        --path_out ${meta.id}_compartments.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
