include { IMAGING_PREPROCESS } from "../../../modules/sanger/imaging/preprocess"
include { IMAGING_GENERATECOMPANION } from "../../../modules/sanger/imaging/generatecompanion"


workflow STITCH {

    take:
    experiments_ch // channel: [ val(meta), [ imaging_experiment ] ]
    psf_folder

    main:
    ch_versions = Channel.empty()

    IMAGING_GENERATECOMPANION(experiments_ch)
    ch_versions = ch_versions.mix(IMAGING_GENERATECOMPANION.out.versions.first())

    tiles = IMAGING_GENERATECOMPANION.out.csv.splitCsv(header: true, strip:true, sep: ",").map {row ->
        [row[0], row[1].index, file(row[1].root_xml), row[1].image_id]
    }
    IMAGING_PREPROCESS(tiles, psf_folder)
    ch_versions = ch_versions.mix(IMAGING_PREPROCESS.out.versions.first())

    emit:
    tiles = IMAGING_PREPROCESS.out.fovs
    companion = IMAGING_GENERATECOMPANION.out.companion

    versions = ch_versions                     // channel: [ versions.yml ]
}
