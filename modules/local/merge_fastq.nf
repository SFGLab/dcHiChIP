
process MERGE_FASTQ {
    tag "$meta.id"

    input:
    tuple val(meta), path(fastq1, stageAs: "?/*"), path(fastq2, stageAs: "?/*")

    output:
    val(meta.id), emit: sample
    val(meta.chipseq), emit: chipseq
    path("${meta.id}_sample_R1.fastq"), emit: fastq1
    path("${meta.id}_sample_R2.fastq"), emit: fastq2
 
    script:
    """
    if test ${fastq1[0]} = 1
    then
        cat 1/${fastq1[1]} > ${meta.id}_sample_R1.fastq
        cat 1/${fastq2[1]} > ${meta.id}_sample_R2.fastq
    else
        cat ${fastq1.join(' ')} > ${meta.id}_sample_R1.fastq
        cat ${fastq2.join(' ')} > ${meta.id}_sample_R2.fastq
    fi
    """
}
