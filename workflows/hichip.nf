/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHichip.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_SE_R1 } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_SE_R2 } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_VIEW as FILTER_QUALITY} from '../modules/nf-core/samtools/view/main'
include { REMOVE_DUPLICATES} from '../modules/local/remove_duplicates'
include { DEEPTOOLS_BAMCOVERAGE } from '../modules/nf-core/deeptools/bamcoverage/main'
include { MACS3_CALLPEAK } from '../modules/nf-core/macs3/callpeak/main'
include { MAPS } from '../modules/local/maps'
include { JUICERTOOLS } from '../modules/local/juicertools'
include { HOMER_ANNOTATEPEAKS } from '../modules/nf-core/homer/annotatepeaks/main'
include { HOMER_FINDMOTIFSGENOME } from '../modules/local/homer/findmotifsgenome/main'
include { DEEPTOOLS_PLOTPCA } from '../modules/nf-core/deeptools/plotpca/main'
include { DEEPTOOLS_PLOTCORRELATION } from '../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTCOVERAGE } from '../modules/local/deeptools/plotcoverage/main'
include { BIOFRAME } from '../modules/local/bioframe'
include { GSTRIPE } from '../modules/local/gstripe'
include { DEEPTOOLS_MUTIBIGWIGSUMMARY } from '../modules/local/deeptools/multibigwigsummary/main'
include { PAIRTOOLS_PARSE2 } from '../modules/local/pairtools/parse2/main'
include { COOLER_CLOAD } from '../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY } from '../modules/nf-core/cooler/zoomify/main'
include { CALDER } from '../modules/local/calder2'
include { COOLTOOLS_INSULATION } from '../modules/local/cooltools/insulation.nf'
include { COOLTOOLS_EIGSCIS } from '../modules/local/cooltools/eigscis.nf'
include { MULTIMM } from '../modules/local/multimm.nf'
include { AWK } from '../modules/local/awk.nf'
include { BEDTOOLS_NUC } from '../modules/nf-core/bedtools/nuc/main'
include { COOLTOOLS_BED_INVERT } from '../modules/local/cooltools_bed_invert.nf'
include { SAMTOOLS_MERGE } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FIXMATE } from '../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_MARKDUP } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_2 } from '../modules/nf-core/samtools/sort/main'
include { FILTER_PAIRES } from '../modules/local/filter_pairs.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow HICHIP {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    Channel
        .fromPath(params.fasta, checkIfExists:true)
        .map{[["id": it.baseName], it]}
        .set{ch_fasta}
    
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
        .fromPath("${params.fasta}.fai", checkIfExists:true)
        .map{[["id": it.baseName], it]}
        .set{ch_fasta_fai}

    Channel
        .fromPath("${params.chrom_size}", checkIfExists:true)
        .set{ch_chrom_size}

    //params.bwa_index ? params.bwa_index : 
    Channel
        .fromPath( params.fasta + '.{bwt,sa,ann,amb,pac}', checkIfExists:true)
        .toSortedList()
        .map{[["id": it[0].baseName], it]}
        .set{ch_bwa_index}
    
    //ch_fasta.view()
    //ch_fasta_fai.view()
    //ch_bwa_index.view()
    def samples_cnt = 0;
    Channel
        .fromPath(params.input, checkIfExists:true)
        .splitCsv(header:true, strip:true)
        .map {samples_cnt += 1; [it, samples_cnt]}
        .multiMap {row, idx -> 
            ch_hichip: tuple(
                ["id": "${row.id}", "single_end": false, "group": row.group, "type": "hichip"], 
                [file(row.hichip_r1, checkIfExists: true),
                file(row.hichip_r2, checkIfExists: true)]
            )
            ch_chipseq: row.chipseq_r1 ? 
                tuple(
                    ["id": "${row.id}", "single_end": false, "group": row.group, "type": "chipseq"], 
                    [file(row.chipseq_r1, checkIfExists: true),
                    file(row.chipseq_r2, checkIfExists: true)]
                ) : tuple()
            ch_narrowpeak: row.narrowpeak ? 
                tuple(
                    ["id": "${row.id}", "single_end": false, "group": row.group], 
                    file(row.narrowpeak, checkIfExists: true)
                ) : tuple()
        }
        .set{ch_input_data}
    
    ch_input_data
    .ch_chipseq
    .filter{it.size() > 0}
    .map{[["id": "${it[0].group}", "single_end": it[0].single_end, "group": it[0].group], it[1]]}
    .groupTuple()
    .map{meta, fastqs -> meta["type"] = "chipseq"; [meta, fastqs.flatten()]}
    .set{ch_chipseq}
    
    ch_input_data
    .ch_hichip
    .map{[["id": "${it[0].group}", "single_end": it[0].single_end, "group": it[0].group], it[1]]}
    .groupTuple()
    .map{meta, fastqs -> meta["type"] = "hichip"; [meta, fastqs.flatten()]}
    .set{ch_hichip}

    ch_input_data
    .ch_narrowpeak
    .filter{it.size() > 0}
    .set{ch_narrowpeak}

    //ch_hichip.view()
    
    CAT_FASTQ(
        ch_hichip
        .filter{(it[0].single_end && it[1].size() > 1) || (!it[0].single_end && it[1].size() > 2)}
        .mix(
            ch_chipseq
            .filter{(it[0].single_end && it[1].size() > 1) || (!it[0].single_end && it[1].size() > 2)}
        )
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    ch_input_data
    .ch_hichip
    .mix(ch_input_data
        .ch_chipseq
        .filter{it.size() > 0}
    )
    .mix(CAT_FASTQ.out.reads)
    .set{ch_bwa_mem_in}
    
    //ch_input_data
    //.ch_hichip.view()
    //ch_bwa_mem_in.view()

    BWA_MEM_SE_R1(
        ch_bwa_mem_in
        .filter{it[0].type == "hichip"}
        .map{
            [it[0], it[1][0]]
        },
        ch_fasta.first(),
        ch_bwa_index,
        true
    )
    BWA_MEM_SE_R2(
        ch_bwa_mem_in
        .filter{it[0].type == "hichip"}
        .map{
            [it[0], it[1][1]]
        },
        ch_fasta.first(),
        ch_bwa_index,
        true
    )
    SAMTOOLS_MERGE(
        BWA_MEM_SE_R1.out.bam
        .join(BWA_MEM_SE_R2.out.bam)
        .map{[it[0], [it[1], it[2]]]},
        ch_fasta.first(),
        ch_fasta_fai.first()
    )
    SAMTOOLS_MERGE.out.bam.view()
    channel.fromPath(["$projectDir/assets/feather", "$projectDir/assets/MAPS"]).toSortedList().view()
    FILTER_PAIRES(
        SAMTOOLS_MERGE.out.bam,
        channel.fromPath(["$projectDir/assets/feather", "$projectDir/assets/MAPS"]).toSortedList()
    )
    SAMTOOLS_FIXMATE(
        FILTER_PAIRES.out.bam
    )
    SAMTOOLS_SORT(
        SAMTOOLS_FIXMATE.out.bam,
        ch_fasta.first()
    )

    SAMTOOLS_MARKDUP(
        SAMTOOLS_SORT.out.bam,
        ch_fasta.first()
    )
    SAMTOOLS_SORT_2(
        SAMTOOLS_MARKDUP.out.bam,
        ch_fasta.first()
    )

    BWA_MEM(
        ch_bwa_mem_in,
        ch_fasta.first(),
        ch_bwa_index,
        false
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    FILTER_QUALITY(
        BWA_MEM.out.bam.map{[it[0], it[1], file("$projectDir/assets/NO_FILE")]},
        ch_fasta.first(),
        []
    )
    ch_versions = ch_versions.mix(FILTER_QUALITY.out.versions)
    
    REMOVE_DUPLICATES(
        FILTER_QUALITY.out.bam,
        ch_fasta.first()
    )
    ch_versions = ch_versions.mix(REMOVE_DUPLICATES.out.versions)

    DEEPTOOLS_BAMCOVERAGE(
        REMOVE_DUPLICATES.out.bam.join(REMOVE_DUPLICATES.out.csi),
        ch_fasta.map{it[1]}.first(),
        ch_fasta_fai.map{it[1]}.first()
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions)

    //REMOVE_DUPLICATES.out.bam.view()
    
    REMOVE_DUPLICATES.out.bam
    .filter{it[0].type == "hichip"}
    .map { [it[0].id, [it[0], it[1]]]}
    .join(
        REMOVE_DUPLICATES.out.bam
        .filter{it[0].type == "chipseq"}
        .map { [it[0].id, [it[0], it[1]]]}
        , remainder: true
    ).join(
        ch_narrowpeak.map{[it[0].id, [it[0], it[1]]]}
        , remainder: true
    ).set{ch_bams}

    MACS3_CALLPEAK(
        ch_bams
        .filter{!it[2] && !it[3]}.map{[it[1][0], it[1][1], []]}
        .mix(
            ch_bams
            .filter{it[2] && !it[3]}.map{[it[2][0], it[2][1], []]}
        ),
        Channel.value(params.genome_size)
    )
    ch_versions = ch_versions.mix(MACS3_CALLPEAK.out.versions)
   
    //ch_bwa_mem_in.view()
    MACS3_CALLPEAK.out.peak
    .map{[it[0].id, it]}
    .mix(ch_narrowpeak.map{[it[0].id, [it[0], it[1]]]})
    .join(ch_bwa_mem_in.filter{it[0].type == "hichip"}
    .map{[it[0].id, it]})
    .join(
        REMOVE_DUPLICATES.out.bam
        .filter{it[0].type == "hichip"}
        .map { [it[0].id, [it[0], it[1]]]}
    )
    .set{ch_maps_fasta_in}
    

    MACS3_CALLPEAK.out.peak
    .map{[it[0].id, it]}
    .mix(ch_narrowpeak.map{[it[0].id, [it[0], it[1]]]})
    .join(
        SAMTOOLS_SORT_2.out.bam
        .map{[it[0].id, it]}
    )
    .set{ch_maps_in}
    

    //ch_maps_fasta_in.view()

    if (!params.skip_maps){
        MAPS(
            ch_maps_fasta_in.map{[it[2][0], it[2][1][0], it[2][1][1]]},
            ch_maps_in.map{[it[2][0], it[2][1]]},
            ch_maps_in.map{it[1]},
            ch_bwa_index,
            ch_genomics_features.first(),
            channel.fromPath(["$projectDir/assets/feather", "$projectDir/assets/MAPS"]).toSortedList()
        )
        ch_versions = ch_versions.mix(MAPS.out.versions)

        JUICERTOOLS(
            MAPS.out.hic
        )
        ch_versions = ch_versions.mix(JUICERTOOLS.out.versions)
        
        GSTRIPE(
            MAPS.out.bedpe
        )
        CALDER(
            JUICERTOOLS.out.hic,
            params.calder_bin,
            params.ref_short
        )

    }
    
    HOMER_ANNOTATEPEAKS(
        ch_maps_in.map{it[1]},
        ch_fasta.map{it[1]}.first(),
        ch_gtf.first()
    )

    HOMER_FINDMOTIFSGENOME(
        ch_maps_in.map{it[1]},
        ch_fasta.map{it[1]}.first()
    )

    BIOFRAME(
        ch_maps_in.map{it[1]},
        ch_jaspar_motif.first(),
        ch_blacklist.first()
    )
    DEEPTOOLS_BAMCOVERAGE.out.bigwig
        .map{[["id": it[0].group, "type":it[0].type], it[0].id, it[1]]}
        .groupTuple()
        .filter{it[2].size() > 1}
        .set{ch_bigwigs}
        
    //ch_bigwigs.view()

    DEEPTOOLS_MUTIBIGWIGSUMMARY(
       ch_bigwigs.map{[it[0], it[2]]},
       ch_bigwigs.map{it[1]}
    )

    DEEPTOOLS_PLOTPCA(
        DEEPTOOLS_MUTIBIGWIGSUMMARY.out.npz
    )
    
    REMOVE_DUPLICATES.out.bam
    .join(REMOVE_DUPLICATES.out.csi)
         .map{[["type":it[0].type, "id": it[0].group], it[0].id, it[1], it[2]]}
         .groupTuple()
         .filter{it[2].size() > 1}
         .set{ch_coverage_in}

    DEEPTOOLS_PLOTCOVERAGE(
        ch_coverage_in.map{[it[0], it[2], it[3]]},
        ch_coverage_in.map{it[1]}
         
    )

    DEEPTOOLS_PLOTCORRELATION(
        DEEPTOOLS_MUTIBIGWIGSUMMARY.out.npz,
        channel.value(params.plot_method),
        channel.value(params.plot_type)
    )

    PAIRTOOLS_PARSE2(
        BWA_MEM.out.bam,
        ch_chrom_size.first()
    )
    
    COOLER_CLOAD(
        PAIRTOOLS_PARSE2.out.nodups.map{[it[0], it[1], [], params.cool_bin]},
        ch_chrom_size.first()

    )

    COOLER_ZOOMIFY(
        COOLER_CLOAD.out.cool.map{[it[0], it[1]]}
    )

    COOLTOOLS_INSULATION(
        COOLER_ZOOMIFY.out.mcool.map{[it[0], it[1]]},
        params.insulation_resultions[params.cooler_zoomify_res]
    )

    COOLTOOLS_EIGSCIS(
        COOLER_ZOOMIFY.out.mcool,
        params.cooler_eigscis_resultion
    )
    
    AWK(
        COOLTOOLS_EIGSCIS.out.scores
    )
    BEDTOOLS_NUC(
        AWK.out.bed.combine(ch_fasta.map{it[1]}).map{[it[0], it[2], it[1]]}
    )
    COOLTOOLS_BED_INVERT(
        AWK.out.bed
    )
    MULTIMM(
        MAPS.out.bedpe,
        COOLTOOLS_BED_INVERT.out.compartments
    )
    
    //
    // MODULE: Run FastQC
    //
    
    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    //
    // MODULE: MultiQC
    //
    /*
    workflow_summary    = WorkflowHichip.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowHichip.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    */
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
