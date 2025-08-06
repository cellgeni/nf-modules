include { IMAGING_PREPROCESS } from "../../../modules/sanger-cellgeni/imaging/preprocess"
include { IMAGING_GENERATECOMPANION } from "../../../modules/sanger-cellgeni/imaging/generatecompanion"

workflow PREPROCESS_TILES {
    take:
    master_file_ch // channel: [ val(meta), val(round_index), file(master_file_folder), val(master_file_name) ]
    psf_folder

    main:
    ch_versions = Channel.empty()

    IMAGING_GENERATECOMPANION(master_file_ch)
    ch_versions = ch_versions.mix(IMAGING_GENERATECOMPANION.out.versions.first())

    tiles = IMAGING_GENERATECOMPANION.out.csv
        .splitCsv(header: true, strip: true, sep: ",")
        .map { row ->
            [row[0], file(row[1].root_xml), row[1].master_file, row[1].image_id, row[1].index]
        }
    IMAGING_PREPROCESS(tiles, psf_folder)
    ch_versions = ch_versions.mix(IMAGING_PREPROCESS.out.versions.first())

    emit:
    companion_tiles = IMAGING_GENERATECOMPANION.out.companion.combine(IMAGING_PREPROCESS.out.fovs, by: 0).groupTuple(by: [0, 1]) // [ meta, companion, [fovs] ]
    versions = ch_versions // channel: [ versions.yml ]
}
