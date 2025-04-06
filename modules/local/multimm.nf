
process MULTIMM {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_multimm:8c581e45aa852599':
        'community.wave.seqera.io/library/pip_multimm:a503e99c16691d6a' }"
    
    input:
    tuple val(meta), path(loops)
    tuple val(meta2), path(peaks)
    tuple val(meta3), path(bwa_index)
    path (force_field)
    path (reqiured_scripts)

    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    tuple val(meta), path("*.hic.input"), emit: hic
    path "versions.yml"                 , emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    
    """
    cat <<-END_CONFIG > config.ini
    [Main]

    ; Platform selection
    PLATFORM = OpenCL

    ; Input data
    FORCEFIELD_PATH = ${force_field}
    LOOPS_PATH = ${loops}
    COMPARTMENT_PATH = /home/skorsak/Data/Rao/subcompartments_primary_replicate/sub_compartments/all_sub_compartments.bed
    ATACSEQ_PATH = /home/skorsak/Data/encode/ATAC-Seq/ENCSR637XSC_GM12878/ENCFF667MDI_pval.bigWig
    OUT_PATH = application_note

    ; Simulation Parameters
    N_BEADS = 50000
    SHUFFLE_CHROMS = True
    NUC_DO_INTERPOLATION = True

    ; Enable forcefield for genome-wide simulation
    SC_USE_SPHERICAL_CONTAINER = True
    CHB_USE_CHROMOSOMAL_BLOCKS = True
    SCB_USE_SUBCOMPARTMENT_BLOCKS = True
    IBL_USE_B_LAMINA_INTERACTION = True
    CF_USE_CENTRAL_FORCE = True

    ; Simulation Parameters
    SIM_RUN_MD = True
    SIM_SAMPLING_STEP = 50
    SIM_N_STEPS = 1000
    TRJ_FRAMES = 100
    END_CONFIG

    MultiMM -c config.ini ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MultiMM: 1.1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.bedpe
    touch ${prefix}.hic.input
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MultiMM: 1.1.0
    END_VERSIONS
    """
}
