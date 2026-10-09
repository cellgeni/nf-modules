#!/usr/bin/env nextflow
//
// adopted from https://github.com/cellgeni/spatialxe/blob/dev/subworkflows/local/proseg_preset_proseg2baysor/main.nf
// Runs proseg for the xenium format and proseg2baysor to generate cell ploygons
//

include { PROSEG                    } from '../../../modules/sanger-cellgeni/proseg/proseg/main'
include { PROSEG_TO_BAYSOR          } from '../../../modules/sanger-cellgeni/proseg/proseg_to_baysor/main'
include { GEOJSON_FEATURECOLLECTION } from '../../../modules/sanger-cellgeni/geojson/featurecollection/main'

workflow PROSEG_PRESET_PROSEG2BAYSOR {
    take:
    _ch_bundle_path // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]

    main:

    ch_versions = channel.empty()
    ch_coordinate_space = channel.value("microns")

    // run proseg with the xenium format
    PROSEG(
        ch_transcripts_parquet,
        "xenium",
        ["parquet", "csv", "mtx.gz"],
    )
    ch_versions = ch_versions.mix(PROSEG.out.versions)

    // run proseg-to-baysor on the data generated with the proseg run
    // PROSEG_TO_BAYSOR(PROSEG.out.transcript_metadata, PROSEG.out.cell_polygons)
    PROSEG_TO_BAYSOR(PROSEG.out.sd_zarr)
    ch_versions = ch_versions.mix(PROSEG_TO_BAYSOR.out.versions)

    GEOJSON_FEATURECOLLECTION(PROSEG_TO_BAYSOR.out.cell_polygons)
    ch_versions = ch_versions.mix(GEOJSON_FEATURECOLLECTION.out.versions)

    emit:
    proseg_sd_zarr   = PROSEG.out.sd_zarr // channel: [ val(meta), [ "cell-polygons.geojson.gz" ] ]
    xr_polygons      = GEOJSON_FEATURECOLLECTION.out.geojson // channel: [ val(meta), [ "xr-cell-polygons.geojson" ] ]
    xr_metadata      = PROSEG_TO_BAYSOR.out.transcript_metadata // channel: [ [ "xr-transcript-metadata.csv" ] ]
    coordinate_space = ch_coordinate_space // channel: [ "microns" ]
    versions         = ch_versions // channel: [ versions.yml ]
}
