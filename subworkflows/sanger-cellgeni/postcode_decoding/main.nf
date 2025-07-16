include { IMAGING_POSTCODE } from '../../../modules/sanger-cellgeni/imaging/postcode/main'
include { IMAGING_POSTCODEPREP } from '../../../modules/sanger-cellgeni/imaging/postcodeprep/main'

workflow POSTCODE_DECODING {
    take:
    ch_profile // channel: [ val(meta), [ profile ] ]
    ch_tabular_codebook // channel: [ val(meta), [ tabular_codebook ] ]
    ch_readout_file // channel: [ val(meta), [ readout_file ] ]
    ch_R // channel: [ val(meta), [ R ] ]
    ch_loc // channel: [ val(meta), [ spot_loc ] ]

    main:

    ch_versions = Channel.empty()

    for_decoding = ch_profile
        .combine(
            ch_tabular_codebook,
            by: 0
        )
        .combine(
            ch_readout_file,
            by: 0
        )
        .combine(ch_R, by: 0)

    for_decoding.view()

    IMAGING_POSTCODEPREP(for_decoding)
    ch_versions = ch_versions.mix(IMAGING_POSTCODEPREP.out.versions.first())

    IMAGING_POSTCODE(IMAGING_POSTCODEPREP.out.for_decoding.combine(ch_loc, by: 0))
    ch_versions = ch_versions.mix(IMAGING_POSTCODEPREP.out.versions.first())

    emit:
    decoding_result = IMAGING_POSTCODE.out.decoded_peaks // channel: [ val(meta), [ pixel/cell ] ]
    processed_profiles = IMAGING_POSTCODEPREP.out.for_decoding // channel: [ val(meta), [ pixel/cell ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
