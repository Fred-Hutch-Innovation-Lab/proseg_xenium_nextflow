/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { PROSEG_XENIUM_ANALYSIS } from '../subworkflows/local/proseg_xenium_analysis'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROSEG_XENIUM_NEXTFLOW {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Prepare input channels from samplesheet
    //
    ch_samplesheet
        .map { meta, xenium_outdir ->
            // Create channels for xenium bundle and transcripts file
            def transcripts_file = file("${xenium_outdir}/transcripts.parquet")
            if (!transcripts_file.exists()) {
                error("Transcripts file not found: ${transcripts_file}")
            }
            return [
                [meta, xenium_outdir],           // xenium bundle channel
                [meta, transcripts_file]         // transcripts channel
            ]
        }
        .transpose()
        .branch { item ->
            xenium_bundle: item[0][1].toString().endsWith(item[0][0].id.toString()) || item[0][1].isDirectory()
                return item
            transcripts: item[1][1].toString().endsWith('transcripts.parquet')
                return item
        }
        .set { ch_inputs }

    // Separate the channels
    ch_xenium_bundle = ch_inputs.xenium_bundle
    ch_transcripts = ch_inputs.transcripts

    //
    // SUBWORKFLOW: Run PROSEG analysis workflow
    //
    PROSEG_XENIUM_ANALYSIS (
        ch_xenium_bundle,
        ch_transcripts,
        'xenium',                                    // proseg mode
        tuple('parquet', 'parquet', 'parquet'),     // output formats
        5,                                          // expansion distance
        [],                                         // coordinate transform (empty)
        [],                                         // nuclei (empty)
        [],                                         // cells (empty)
    )
    ch_versions = ch_versions.mix(PROSEG_XENIUM_ANALYSIS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'proseg_xenium_nextflow_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
