include { IMAGINGGENERATECOMPANION } from "../modules/sanger/imaging/generatecompanion"
include { IMAGINGPREPROCESS } from "../modules/sanger/imaging/preprocess"


workflow STITCH {

    take:
    experiments_ch // channel: [ val(meta), [ imaging_experiment ] ]
    psf_folder

    main:
    ch_versions = Channel.empty()

    IMAGINGGENERATECOMPANION(experiments_ch)
    ch_versions = ch_versions.mix(IMAGINGGENERATECOMPANION.out.versions.first())

    tiles = IMAGINGGENERATECOMPANION.out.csv.splitCsv(header: true, strip:true, sep: ",").map {row ->
        [row[0], row[1].index, file(row[1].root_xml), row[1].image_id]
    }
    IMAGINGPREPROCESS(tiles, psf_folder)
    ch_versions = ch_versions.mix(IMAGINGPREPROCESS.out.versions.first())

    emit:
    tiles = IMAGINGPREPROCESS.out.fovs
    companion = IMAGINGGENERATECOMPANION.out.companion

    versions = ch_versions                     // channel: [ versions.yml ]
}