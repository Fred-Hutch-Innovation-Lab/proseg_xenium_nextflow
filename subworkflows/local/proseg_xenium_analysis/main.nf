//
// Subworkflow to run PROSEG, PROSEG_TO_BAYSOR, and XENIUMRANGER_IMPORT_SEGMENTATION
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROSEG                           } from '../../../modules/nf-core/proseg/proseg/main'
include { PROSEG_TO_BAYSOR                 } from '../../../modules/nf-core/proseg/proseg_to_baysor/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROSEG_XENIUM_ANALYSIS {

    take:
    ch_xenium_bundle                    // channel: [ val(meta), path(xenium_bundle) ]
    ch_transcripts                      // channel: [ val(meta), path(transcripts) ]
    proseg_mode                         // val: proseg mode (e.g., 'xenium')
    output_formats                      // tuple: (transcript_metadata_fmt, cell_metadata_fmt, expected_counts_fmt)
    expansion_distance                  // val: expansion distance for xeniumranger
    coordinate_transform                // path: coordinate transform file (optional)
    nuclei                             // path: nuclei file (optional)
    cells                              // path: cells file (optional)

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run PROSEG to generate cell segmentation
    //
    PROSEG (
        ch_transcripts,
        proseg_mode,
        output_formats
    )
    ch_versions = ch_versions.mix(PROSEG.out.versions)

    //
    // MODULE: Convert PROSEG output to Baysor format
    //
    PROSEG_TO_BAYSOR (
        PROSEG.out.transcript_metadata,
        PROSEG.out.cell_polygons
    )
    ch_versions = ch_versions.mix(PROSEG_TO_BAYSOR.out.versions)

    //
    // MODULE: Import segmentation into Xenium Ranger
    //
    XENIUMRANGER_IMPORT_SEGMENTATION (
        ch_xenium_bundle,
        expansion_distance,
        coordinate_transform,
        nuclei,
        cells,
        PROSEG_TO_BAYSOR.out.baysor_transcript_metadata,
        PROSEG_TO_BAYSOR.out.baysor_cell_polygons
    )
    ch_versions = ch_versions.mix(XENIUMRANGER_IMPORT_SEGMENTATION.out.versions)

    emit:
    // PROSEG outputs
    proseg_transcript_metadata     = PROSEG.out.transcript_metadata
    proseg_cell_polygons          = PROSEG.out.cell_polygons
    proseg_cell_metadata          = PROSEG.out.cell_metadata
    proseg_cell_polygons_layers   = PROSEG.out.cell_polygons_layers
    proseg_expected_counts        = PROSEG.out.expected_counts
    proseg_union_cell_polygons    = PROSEG.out.union_cell_polygons
    proseg_maxpost_counts         = PROSEG.out.maxpost_counts
    proseg_output_rates           = PROSEG.out.output_rates
    proseg_cell_hulls             = PROSEG.out.cell_hulls
    proseg_gene_metadata          = PROSEG.out.gene_metadata
    proseg_metagene_rates         = PROSEG.out.metagene_rates
    proseg_metagene_loadings      = PROSEG.out.metagene_loadings
    proseg_cell_voxels            = PROSEG.out.cell_voxels

    // PROSEG_TO_BAYSOR outputs
    baysor_cell_polygons          = PROSEG_TO_BAYSOR.out.baysor_cell_polygons
    baysor_transcript_metadata    = PROSEG_TO_BAYSOR.out.baysor_transcript_metadata

    // XENIUMRANGER outputs
    xeniumranger_outs             = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs

    // Versions
    versions                      = ch_versions
}