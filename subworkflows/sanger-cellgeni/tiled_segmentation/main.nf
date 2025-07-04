include { IMAGING_CELLPOSE as CELLPOSE } from '../../../modules/sanger/imaging/cellpose/main'
include { IMAGING_STARDIST as STARDIST } from '../../../modules/sanger/imaging/stardist/main'
include { IMAGING_INSTANSEG as INSTANSEG } from '../../../modules/sanger/imaging/instanseg/main'
include { IMAGING_DEEPCELL as DEEPCELL } from '../../../modules/sanger/imaging/deepcell/main'
include { IMAGING_MERGEOUTLINES as MERGEOUTLINES } from '../../../modules/sanger/imaging/mergeoutlines/main'
include { IMAGING_GENERATETILES as GENERATE_TILE_COORDS } from '../../../modules/sanger/imaging/generatetiles/main'


workflow TILED_SEGMENTATION {
    take:
    images
    method

    main:
    ch_versions = Channel.empty()
    GENERATE_TILE_COORDS(images)
    ch_versions = ch_versions.mix(GENERATE_TILE_COORDS.out.versions.first())

    images_tiles = GENERATE_TILE_COORDS.out.tile_coords
        .splitCsv(header: true, sep: ",")
        .map { meta, coords ->
            [meta, coords.X_MIN, coords.Y_MIN, coords.X_MAX, coords.Y_MAX]
        }
    tiles_and_images = images_tiles.combine(images, by: 0)
    if (method == "CELLPOSE") {
        CELLPOSE(tiles_and_images.combine(channel.from(params.cell_diameters)))
        wkts = CELLPOSE.out.wkts.groupTuple(by: 0)
        ch_versions = ch_versions.mix(CELLPOSE.out.versions.first())
    }
    else if (method == "STARDIST") {
        STARDIST(tiles_and_images)
        wkts = STARDIST.out.wkts.groupTuple(by: 0)
        ch_versions = ch_versions.mix(STARDIST.out.versions.first())
    }
    else if (method == "INSTANSEG") {
        INSTANSEG(tiles_and_images)
        wkts = INSTANSEG.out.wkts.groupTuple(by: 0)
        ch_versions = ch_versions.mix(INSTANSEG.out.versions.first())
    }
    else if (method == "DEEPCELL") {
        DEEPCELL(tiles_and_images)
        wkts = DEEPCELL.out.wkts.groupTuple(by: 0)
        ch_versions = ch_versions.mix(DEEPCELL.out.versions.first())
    }
    else {
        error("Invalid segmentation method: ${method}. Expected one of: CELLPOSE, STARDIST, INSTANSEG, DEEPCELL")
    }
    MERGEOUTLINES(wkts)
    ch_versions = ch_versions.mix(MERGEOUTLINES.out.versions.first())

    emit:
    geojson = MERGEOUTLINES.out.multipoly_geojsons // channel: [ val(meta), [ geojson ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
