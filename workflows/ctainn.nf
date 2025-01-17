/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                              } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ASTAIR                        	                            } from '../modules/local/astair/main'
include { FASTQC                                                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                                       } from '../modules/nf-core/multiqc/main'
include { PYMBIAS                                                       } from '../modules/local/pymbias/main'
include { RASTAIR                       	                            } from '../modules/local/rastair/main'
include { BWAMEM2                                                       } from '../modules/local/bwamem2/main'
include { BWA_METH                 	                                    } from '../modules/local/bwameth/main'
include { CAT_FASTQ                     	                            } from '../modules/nf-core/cat/fastq/main'
include { RASTAIR_MBIAS                 	                            } from '../modules/local/rastair/main'
include { SAMTOOLS_INDEX                	                            } from '../modules/nf-core/samtools/index/main'
include { RASTAIR_PERREAD               	                            } from '../modules/local/rastair/main'
include { RASTAIR_SUMMARY               	                            } from '../modules/local/rastair/main'
include { GATK4_MARKDUPLICATES          	                            } from '../modules/local/gatk4/markduplicates/main'
include { PREP_BEDTOOLS_INTERSECT       	                            } from '../modules/local/bedtools/prep_bedtools_intersect' 
include { ASTAIR_BEDTOOLS_INTERSECT     	                            } from '../modules/local/bedtools/astair_bedtools_intersect' 
include { RASTAIR_BEDTOOLS_INTERSECT    	                            } from '../modules/local/bedtools/rastair_bedtools_intersect' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                        } from '../subworkflows/local/utils_nfcore_ctainn_pipeline'
include { validateInputSamplesheet                                      } from '../subworkflows/local/utils_nfcore_ctainn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                            CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_metbed                                                               = Channel.fromPath(params.metbed).map { it -> [[id:it.Name], it] }.collect()
ch_metdir                                                               = Channel.fromPath(params.metdir).map { it -> [[id:it.Name], it] }.collect()
ch_metref                                                               = Channel.fromPath(params.metref).map { it -> [[id:it.Name], it] }.collect()
ch_metfai                                                               = Channel.fromPath(params.metfai).map { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                               IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

def getMbiasParams(tsv_file) {
    def values = [:]
    def content = file(tsv_file).readLines()
    def lines = content[1..-1]
    def ot_values = []
    def ob_values = []
    
    lines.each { line ->
        def fields = line.split('\t')
        def strand = fields[1]
        def left = fields[3]
        def right = fields[4]
        
        if (strand == "OT") {
            ot_values.add(left)
            ot_values.add(right)
        } else if (strand == "OB") {
            ob_values.add(left)
            ob_values.add(right)
        }
    }
    
    values.not = ot_values.join(',')
    values.nob = ob_values.join(',')
    
    return values
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                              RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CTAINN {

    take:
    ch_samplesheet
    
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Create branch for single and/or multiple lanes.
    //
    ch_samplesheet
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }
    ch_fastq_single     = ch_fastq.single
    ch_fastq_multiple   = ch_fastq.multiple

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (ch_cat_fastq)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    // User aligner selection
    if (params.aligner == 'bwa-meth') {
        //
        // MODULE: Run BWA-METH
        //
        BWA_METH(ch_cat_fastq, ch_metdir, ch_metref)
        ch_versions = ch_versions.mix(BWA_METH.out.versions.first())
        ch_bam_mapped = BWA_METH.out.bam
    } else if (params.aligner == 'bwa-mem2') {
        //
        // MODULE: Run BWA-MEM2
        //
        sort_bam = 'sort'
        BWAMEM2(ch_cat_fastq, ch_metdir, ch_metref, ch_metfai, sort_bam, sort_bam)
        ch_versions = ch_versions.mix(BWAMEM2.out.versions.first())
        ch_bam_mapped = BWAMEM2.out.bam
    } else {
        error "Invalid aligner selected: ${params.aligner}. Please choose either 'bwa-meth' or 'bwa-mem2'"
    }

    //
    // MODULE: Run GATK4 MarkDuplicates
    //
    GATK4_MARKDUPLICATES(ch_bam_mapped, params.metref, params.metfai)
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
    ch_bam_dedup = GATK4_MARKDUPLICATES.out.bam

    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    PREP_BEDTOOLS_INTERSECT(ch_bam_dedup, ch_metbed, 'targeted')
    ch_versions = ch_versions.mix(PREP_BEDTOOLS_INTERSECT.out.versions.first())
    ch_bam_mapped_targeted = PREP_BEDTOOLS_INTERSECT.out.bam

    //
    // MODULE: Run SamTools Index
    //
    SAMTOOLS_INDEX(ch_bam_mapped_targeted)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_bam_mapped_targeted_indexed = PREP_BEDTOOLS_INTERSECT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    //
    // MODULE: Run rasTair in per-read mode
    //
    RASTAIR_MBIAS(ch_bam_mapped_targeted_indexed, ch_metdir)
    ch_versions = ch_versions.mix(RASTAIR_MBIAS.out.versions.first())

    //
    // MODULE: Run PyMbias
    //
    PYMBIAS(RASTAIR_MBIAS.out.mbias, params.pymbias_plot_type, params.pymbias_plot_ax_x, params.pymbias_plot_ax_y)
    ch_versions = ch_versions.mix(PYMBIAS.out.versions.first())
    PYMBIAS.out.cutoffs
        .map { meta, cutoff_file -> getMbiasParams(cutoff_file) }
        .set { ch_cutoffs }

    //
    // MODULE: Run rasTair
    //
    RASTAIR(ch_bam_mapped_targeted_indexed, ch_metdir, ch_cutoffs)
    ch_versions = ch_versions.mix(RASTAIR.out.versions.first())
    ch_rastair_mods = RASTAIR.out.mods

    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    RASTAIR_BEDTOOLS_INTERSECT(ch_rastair_mods, ch_metbed, '_rastair_output_targeted')
    ch_versions = ch_versions.mix(RASTAIR_BEDTOOLS_INTERSECT.out.versions.first())
    ch_rastair_mods_targeted = RASTAIR_BEDTOOLS_INTERSECT.out.mods

    //
    // MODULE: Run rasTair summary for the targeted mods
    //
    RASTAIR_SUMMARY(ch_rastair_mods_targeted)
    ch_versions = ch_versions.mix(RASTAIR_BEDTOOLS_INTERSECT.out.versions.first())

    //
    // MODULE: Run astair
    //
    ASTAIR(ch_bam_mapped_targeted_indexed, ch_metdir, params.read_length)
    ch_versions = ch_versions.mix(ASTAIR.out.versions.first())
    ch_astair_mods = ASTAIR.out.mods

    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    ASTAIR_BEDTOOLS_INTERSECT(ch_astair_mods, ch_metbed, '_astair_output_targeted')
    ch_versions = ch_versions.mix(ASTAIR_BEDTOOLS_INTERSECT.out.versions.first())
    ch_astair_mods_targeted = ASTAIR_BEDTOOLS_INTERSECT.out.mods

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
