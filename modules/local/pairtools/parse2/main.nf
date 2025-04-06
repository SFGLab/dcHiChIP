process PAIRTOOLS_PARSE2 {
    tag "$meta.id"
    label 'process_high'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0' :
        'biocontainers/pairtools:1.0.2--py39h2a9f597_0' }"

    input:
    tuple val(meta), path(bam)
    path chromsizes

    output:
    tuple val(meta), path("*.nodups.pairs.gz")  , emit: nodups
    tuple val(meta), path("*.dups.pairs.gz"), emit: dups
    tuple val(meta), path("*.unmapped.pairs.gz"), emit: unmapped
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: '3'
    def args3 = task.ext.args3 ?: '8'
    def args4 = task.ext.args4 ?: '3'
    def args5 = task.ext.args5 ?: '3'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pairtools parse2 \\
        $args \\
        --nproc-in $args2 \\
        --nproc-out $args3  \\
        -c $chromsizes \\
        $bam | \\
        pairtools sort --tmpdir ./ --nproc $args4  | \\
        pairtools dedup --n-proc $args5  \\
        --output ${prefix}.nodups.pairs.gz \\
        --output-dups ${prefix}.dups.pairs.gz \\
        --output-unmapped ${prefix}.unmapped.pairs.gz 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.nodups.pairs.gz
    touch ${prefix}.dups.pairs.gz
    touch ${prefix}.unmapped.pairs.gz 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
}
