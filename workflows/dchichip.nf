/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_dchichip_pipeline'


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


workflow DCHICHIP {

    take:
    ch_fasta
    ch_fasta_fai
    ch_bwa_index
    ch_genomics_features
    ch_jaspar_motif
    ch_blacklist
    ch_gtf
    ch_chrom_size
    ch_chipseq
    ch_hichip
    ch_narrowpeak
    ch_genome_size
    skip_maps
    calder_bin
    ref_short
    calder_chrom
    plot_method
    plot_type
    cool_bin
    insulation_resultions
    cooler_eigscis_resultion
    ch_multiqc_config
    ch_multiqc_custom_config
    ch_multiqc_logo
    ch_multiqc_methods_description
    outdir

    main:
    ch_versions = Channel.empty()
    ch_multimm_in = Channel.empty()

    ch_chipseq
    .map{[["id": "${it[0].group}", "single_end": it[0].single_end, "group": it[0].group], it[1]]}
    .groupTuple()
    .map{meta, fastqs -> meta["type"] = "chipseq"; [meta, fastqs.flatten()]}
    .set{ch_chipseq_grouped}

    ch_hichip
    .map{[["id": "${it[0].group}", "single_end": it[0].single_end, "group": it[0].group], it[1]]}
    .groupTuple()
    .map{meta, fastqs -> meta["type"] = "hichip"; [meta, fastqs.flatten()]}
    .set{ch_hichip_grouped}

    FASTQC(
        ch_chipseq
        .mix(ch_hichip)
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    CAT_FASTQ(
        ch_hichip_grouped
        .filter{(it[0].single_end && it[1].size() > 1) || (!it[0].single_end && it[1].size() > 2)}
        .mix(
            ch_chipseq_grouped
            .filter{(it[0].single_end && it[1].size() > 1) || (!it[0].single_end && it[1].size() > 2)}
        )
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    ch_hichip
    .mix(ch_chipseq)
    .mix(CAT_FASTQ.out.reads)
    .set{ch_bwa_mem_in}

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
    ch_versions = ch_versions.mix(BWA_MEM_SE_R1.out.versions)

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
    ch_versions = ch_versions.mix(BWA_MEM_SE_R2.out.versions)
    SAMTOOLS_MERGE(
        BWA_MEM_SE_R1.out.bam
        .join(BWA_MEM_SE_R2.out.bam)
        .map{[it[0], [it[1], it[2]]]},
        ch_fasta.first(),
        ch_fasta_fai.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    FILTER_PAIRES(
        SAMTOOLS_MERGE.out.bam,
        channel.fromPath(["$projectDir/assets/feather", "$projectDir/assets/MAPS"]).toSortedList()
    )
    ch_versions = ch_versions.mix(FILTER_PAIRES.out.versions)

    SAMTOOLS_FIXMATE(
        FILTER_PAIRES.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions)

    SAMTOOLS_SORT(
        SAMTOOLS_FIXMATE.out.bam,
        ch_fasta.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    SAMTOOLS_MARKDUP(
        SAMTOOLS_SORT.out.bam,
        ch_fasta.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    SAMTOOLS_SORT_2(
        SAMTOOLS_MARKDUP.out.bam,
        ch_fasta.first()
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_2.out.versions)

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
        ch_genome_size
    )
    ch_versions = ch_versions.mix(MACS3_CALLPEAK.out.versions)


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


    if (!skip_maps){
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
        ch_versions = ch_versions.mix(GSTRIPE.out.versions)

        CALDER(
            JUICERTOOLS.out.hic,
            calder_bin,
            ref_short,
            calder_chrom
        )
        ch_versions = ch_versions.mix(CALDER.out.versions)
        ch_multimm_in = ch_multimm_in.mix(MAPS.out.bedpe)
    }

    HOMER_ANNOTATEPEAKS(
        ch_maps_in.map{it[1]},
        ch_fasta.map{it[1]}.first(),
        ch_gtf.first()
    )
    ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)

    HOMER_FINDMOTIFSGENOME(
        ch_maps_in.map{it[1]},
        ch_fasta.map{it[1]}.first()
    )
    ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME.out.versions)

    BIOFRAME(
        ch_maps_in.map{it[1]},
        ch_jaspar_motif.first(),
        ch_blacklist.first()
    )
    ch_versions = ch_versions.mix(BIOFRAME.out.versions)

    DEEPTOOLS_BAMCOVERAGE.out.bigwig
        .map{[["id": it[0].group, "type":it[0].type], it[0].id, it[1]]}
        .groupTuple()
        .filter{it[2].size() > 1}
        .set{ch_bigwigs}

    DEEPTOOLS_MUTIBIGWIGSUMMARY(
       ch_bigwigs.map{[it[0], it[2]]},
       ch_bigwigs.map{it[1]}
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_MUTIBIGWIGSUMMARY.out.versions)

    DEEPTOOLS_PLOTPCA(
        DEEPTOOLS_MUTIBIGWIGSUMMARY.out.npz
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPCA.out.versions)

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
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCOVERAGE.out.versions)

    DEEPTOOLS_PLOTCORRELATION(
        DEEPTOOLS_MUTIBIGWIGSUMMARY.out.npz,
        plot_method,
        plot_type
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_PLOTCORRELATION.out.versions)

    PAIRTOOLS_PARSE2(
        BWA_MEM.out.bam.filter{it[0].type == "hichip"},
        ch_chrom_size.first()
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_PARSE2.out.versions)

    COOLER_CLOAD(
        PAIRTOOLS_PARSE2.out.nodups.map{[it[0], it[1], [], cool_bin]},
        ch_chrom_size.first()
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    COOLER_ZOOMIFY(
        COOLER_CLOAD.out.cool.map{[it[0], it[1]]}
    )
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    COOLTOOLS_INSULATION(
        COOLER_ZOOMIFY.out.mcool.map{[it[0], it[1]]},
        insulation_resultions
    )
    ch_versions = ch_versions.mix(COOLTOOLS_INSULATION.out.versions)

    COOLTOOLS_EIGSCIS(
        COOLER_ZOOMIFY.out.mcool,
        cooler_eigscis_resultion
    )
    ch_versions = ch_versions.mix(COOLTOOLS_EIGSCIS.out.versions)

    AWK(
        COOLTOOLS_EIGSCIS.out.scores
    )
    ch_versions = ch_versions.mix(AWK.out.versions)

    BEDTOOLS_NUC(
        AWK.out.bed.combine(ch_fasta.map{it[1]}).map{[it[0], it[2], it[1]]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_NUC.out.versions)

    COOLTOOLS_BED_INVERT(
        BEDTOOLS_NUC.out.bed
    )
    ch_versions = ch_versions.mix(COOLTOOLS_BED_INVERT.out.versions)

    MULTIMM(
        ch_multimm_in,
        COOLTOOLS_BED_INVERT.out.compartments
    )
    ch_versions = ch_versions.mix(MULTIMM.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'proteinfold_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description   = Channel.value(methodsDescriptionText(ch_multiqc_methods_description))

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files
                        .mix(FASTQC.out.zip.map{it[1]})

    MULTIQC (
        ch_multiqc_files.collect(sort: true),
        ch_multiqc_config.collect()
            .ifEmpty([]),
        ch_multiqc_custom_config
            .collect()
            .ifEmpty([]),
        ch_multiqc_logo
            .collect()
            .ifEmpty([]),
        [],
        []
    )

    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    emit:
    versions        = ch_versions
    multiqc_report  = MULTIQC.out.report
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
