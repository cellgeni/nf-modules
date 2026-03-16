include { OMEZARR_PREPROCESS } from "../../../modules/sanger-cellgeni/omezarr/preprocess"
include { IMAGING_GENERATECOMPANION } from "../../../modules/sanger-cellgeni/imaging/generatecompanion"

workflow PREPROCESS_OME_ZARR_TILES {
    take:
    master_file_ch // channel: [ val(meta), file(master_file_folder) ]
    psf_folder

    main:
    ch_versions = Channel.empty()

    IMAGING_GENERATECOMPANION(master_file_ch)
    ch_versions = ch_versions.mix(IMAGING_GENERATECOMPANION.out.versions.first())

    tiles = IMAGING_GENERATECOMPANION.out.csv
        .splitCsv(header: true, strip: true, sep: ",")
        .map { row ->
            [row[0], file(row[1].root_xml), row[1].image_id, row[1].path]
        }
    OMEZARR_PREPROCESS(tiles, psf_folder)
    ch_versions = ch_versions.mix(OMEZARR_PREPROCESS.out.versions.first())

    emit:
    companion_tiles = IMAGING_GENERATECOMPANION.out.companion.combine(OMEZARR_PREPROCESS.out.fovs, by: 0).groupTuple(by: [0, 1]) // [ meta, companion, [fovs] ]
    versions = ch_versions // channel: [ versions.yml ]
}
