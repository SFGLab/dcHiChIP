process FILTER_PAIRES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "community.wave.seqera.io/library/bedtools_bwa_pybedtools_pysam_pruned:8a65a012e7f8bbf1"

    input:
    tuple val(meta), path(bam)
    path (scripts)

    output:
    tuple val(meta), path('qc_file.txt')       , emit: qc
    tuple val(meta), path('*.bam')       , emit: bam
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in sfglab/hichip/bin/
    """
    run_maps_filter_pairs.sh \\
    --outdir ./ \\
    --qc-file qc_file.txt \\
    --input-bam ${bam} \\
    --command filter_pair

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    """
    touch qc_file.txt
    touch ${meta.id}_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
