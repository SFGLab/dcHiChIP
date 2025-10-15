process AWK {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta) , path(egis_scores)

    output:
    tuple val(meta), path("*_temp.bed"),  emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in sfglab/dchichip/bin/
    """
    # converting tsv file to bed format:
    awk 'NR > 1 && \$1 != "" && \$2 != "" && \$3 != "" && \$5 != "" {OFS="\\t"; print \$1, \$2, \$3, \$5}' ${egis_scores} > ${egis_scores.name}_temp.bed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:

    """
    touch ${egis_scores.name}_temp.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
