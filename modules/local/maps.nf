
process MAPS {
    tag "$meta.id"
    label 'process_medium'
    container "community.wave.seqera.io/library/bedtools_bwa_pybedtools_pysam_pruned:8a65a012e7f8bbf1"
    
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
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`    
    run_maps_pipeline.sh \\
    -i ${prefix} \\
    -p ${peaks} \\
    -t ${task.cpus} \\
    -x \$INDEX \\
    -1 ${fastq1} \\
    -2 ${fastq2} \\
    -u ${features} \\
    ${args}
    
    cp MAPS_output/${prefix}_current/${prefix}.5k.2.sig3Dinteractions.bedpe ${prefix}.bedpe
    cp feather_output/${prefix}_current/${prefix}.hic.input .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
