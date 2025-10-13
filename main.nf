#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sfglab/dchichip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sfglab/dchichip
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DCHICHIP } from './workflows/dchichip'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_dchichip_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_dchichip_pipeline'


//
// WORKFLOW: Run main sfglab/dchichip analysis pipeline
//
workflow SFGLAB_DCHICHIP {
    take:
    samplesheet
    main:
    ch_samplesheet              = samplesheet
    ch_multiqc                  = Channel.empty()
    ch_versions                 = Channel.empty()
    multiqc_report = Channel.empty()
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ).first()  : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo ).first()    : Channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    
    //
    // WORKFLOW: Run pipeline
    //

    DCHICHIP (
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_methods_description,
        params.outdir
    )
    emit:
    multiqc_report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )
    //
    // WORKFLOW: Run main workflow
    //
    SFGLAB_DCHICHIP (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        SFGLAB_DCHICHIP.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
