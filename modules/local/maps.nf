
process MAPS {
    tag "$meta.id"
    container "community.wave.seqera.io/library/bedtools_pybedtools_pysam_samtools_pruned:442b9f4bf56bd920"
    
    input:
    tuple val(meta), path(fastq1), path(fastq2)
    tuple val(meta2), path(peaks)
    tuple val(meta3), path(bwa_index)
    path (features)
    path (reqiured_scripts)

    output:
    tuple val(meta), path("${sample}.bedpe"), emit: bedpe
    tuple val(meta), path("${sample}.hic.input"), emit: hic
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    """
}
