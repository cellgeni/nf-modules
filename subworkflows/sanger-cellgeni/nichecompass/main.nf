include { SPATIAL_NICHECOMPASSTRAINING } from '../../../modules/sanger-cellgeni/spatial/nichecompasstraining/main'
include { SPATIAL_NICHECOMPASSANALYSIS } from '../../../modules/sanger-cellgeni/spatial/nichecompassanalysis/main'

workflow NICHECOMPASS {
    take:
    ch_h5ad // channel: [ val(meta), [ h5ad ] ]

    main:
    ch_versions = Channel.empty()

    SPATIAL_NICHECOMPASSTRAINING(ch_h5ad)
    ch_versions = ch_versions.mix(SPATIAL_NICHECOMPASSTRAINING.out.versions.first())

    SPATIAL_NICHECOMPASSANALYSIS(SPATIAL_NICHECOMPASSTRAINING.out.nichecompass_model)
    ch_versions = ch_versions.mix(SPATIAL_NICHECOMPASSANALYSIS.out.versions.first())

    emit:
    nichecompass_dir = SPATIAL_NICHECOMPASSANALYSIS.out.nichecompass_dir // channel: [ val(meta), [ nichecompass_dir ] ]
    notebook = SPATIAL_NICHECOMPASSANALYSIS.out.notebook // channel: [ val(meta), [ notebook.ipynb ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
