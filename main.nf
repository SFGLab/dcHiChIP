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
    multiqc_report              = Channel.empty()
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true).first()
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ).first()  : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo ).first()    : Channel.empty()
    ch_multiqc_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    Channel
        .fromPath(params.fasta, checkIfExists:true)
        .map{[["id": it.baseName], it]}
        .set{ch_fasta}

    Channel
        .fromPath("${params.fasta}.fai", checkIfExists:true)
        .map{[["id": it.baseName], it]}
        .set{ch_fasta_fai}
    Channel
        .fromPath( params.fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true)
        .toSortedList()
        .map{[["id": it[0].baseName], it]}
        .set{ch_bwa_index}

    Channel
        .fromPath(params.genomics_features, checkIfExists:true)
        .set{ch_genomics_features}

    Channel
        .fromPath(params.jaspar_motif, checkIfExists:true)
        .set{ch_jaspar_motif}

    Channel
        .fromPath(params.blacklist, checkIfExists:true)
        .set{ch_blacklist}

    Channel
        .fromPath(params.gtf, checkIfExists:true)
        .set{ch_gtf}

    Channel
        .fromPath("${params.chrom_size}", checkIfExists:true)
        .set{ch_chrom_size}

    ch_samplesheet
        .multiMap { meta, hichip_r1, hichip_r2, chipseq_r1, chipseq_r2, narrowpeak ->
            hichip: tuple(
                ["id": "${meta.id}", "single_end": false, "group": "${meta.group}", "type": "hichip"],
                [hichip_r1, hichip_r2]
            )
            chipseq: chipseq_r1 ?
                tuple(
                    ["id": "${meta.id}", "single_end": false, "group": "${meta.group}", "type": "chipseq"],
                    [chipseq_r1, chipseq_r2]
                ) : tuple()
            narrowpeak: narrowpeak ?
                tuple(
                    ["id": "${meta.id}", "single_end": false, "group": "${meta.group}"],
                    narrowpeak
                ) : tuple()
        }
        .set{ch_input_data}

    //
    // WORKFLOW: Run pipeline
    //

    DCHICHIP (
        ch_fasta,
        ch_fasta_fai,
        ch_bwa_index,
        ch_genomics_features,
        ch_jaspar_motif,
        ch_blacklist,
        ch_gtf,
        ch_chrom_size,
        ch_input_data.chipseq.filter{it.size() > 0},
        ch_input_data.hichip,
        ch_input_data.narrowpeak.filter{it.size() > 0},
        params.genome_size,
        params.skip_maps,
        params.calder_bin,
        params.ref_short,
        params.calder_chrom,
        params.plot_method,
        params.plot_type,
        params.cool_bin,
        params.insulation_resultions[params.cooler_zoomify_res],
        params.cooler_eigscis_resultion,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_methods_description,
        params.outdir
    )

    multiqc_report = multiqc_report.mix(DCHICHIP.out.multiqc_report)

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
