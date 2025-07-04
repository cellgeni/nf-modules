#!/usr/bin/env/ nextflow
include { IMAGING_GENERATETILES as GENERATE_TILE_COORDS } from '../../../modules/sanger/imaging/generatetiles/main'
include { IMAGING_TILEDSPOTIFLOW as SPOTIFLOW           } from '../../../modules/sanger/imaging/tiledspotiflow/main'
include { IMAGING_MERGEPEAKS as MERGEPEAKS              } from '../../../modules/sanger/imaging/mergepeaks/main'
include { IMAGING_CONCATENATEWKTS as CONCATENATEWKTS    } from '../../../modules/sanger/imaging/concatenatewkts/main'


workflow TILED_SPOTIFLOW {
    take:
    images
    chs_to_call_peaks

    main:
    ch_versions = Channel.empty()
    GENERATE_TILE_COORDS(images)
    images_tiles = GENERATE_TILE_COORDS.out.tile_coords
        .splitCsv(header: true, sep: ",")
        .map { meta, coords ->
            [meta, coords.X_MIN, coords.Y_MIN, coords.X_MAX, coords.Y_MAX]
        }

    SPOTIFLOW(images_tiles.combine(images, by: 0).combine(chs_to_call_peaks))
    ch_versions = ch_versions.mix(SPOTIFLOW.out.versions.first())

    MERGEPEAKS(SPOTIFLOW.out.peaks.groupTuple(by: [0, 1]))
    ch_versions = ch_versions.mix(MERGEPEAKS.out.versions.first())

    CONCATENATEWKTS(MERGEPEAKS.out.merged_peaks.groupTuple(by: 0))
    ch_versions = ch_versions.mix(CONCATENATEWKTS.out.versions.first())

    emit:
    spots_csv = CONCATENATEWKTS.out.concatenated_peaks
    versions  = ch_versions
}
