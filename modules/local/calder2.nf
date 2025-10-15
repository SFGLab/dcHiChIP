
process CALDER {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/r-calder2_r-r6:62ef6858209773b9':
        'community.wave.seqera.io/library/r-calder2_r-r6:62ef6858209773b9' }"

    input:
    tuple val(meta), path(hic)
    val(bin_size)
    val(genome_build)
    val(chrom_ls)

    output:
    tuple val(meta), path("*_output"), emit: outputs
    path "versions.yml"                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def version = "0.3"
    """
    #!/usr/bin/env Rscript

    # and load Calder
    library(CALDER)
    chrs = c(${chrom_ls})
    CALDER(contact_file_hic="${hic}",
			chrs=chrs,
			bin_size=${bin_size},
			genome='${genome_build}',
			save_dir='${prefix}_output',
			save_intermediate_data=FALSE,
			n_cores=${task.cpus},
			sub_domains=FALSE)

    writeLines(c('"${task.process}":', '    CALDER: $version'), "versions.yml")

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "0.3"
    """
    mkdir ${prefix}_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CALDER: $version
    END_VERSIONS
    """
}
