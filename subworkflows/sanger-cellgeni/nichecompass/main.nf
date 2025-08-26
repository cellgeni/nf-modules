include { NICHECOMPASS_TRAINING } from '../../../modules/sanger-cellgeni/imaging/nichecompass/training/main'
include { NICHECOMPASS_ANALYSIS } from '../../../modules/sanger-cellgeni/imaging/nichecompass/analysis/main'

workflow NICHECOMPASS {

    take:
    ch_h5ad // channel: [ val(meta), [ h5ad ] ]

    main:

    ch_versions = Channel.empty()

    NICHECOMPASS_TRAINING (ch_h5ad)
    ch_versions = ch_versions.mix(NICHECOMPASS_TRAINING.out.versions.first())

    NICHECOMPASS_ANALYSIS (NICHECOMPASS_TRAINING.out.nichecompass_model)
    ch_versions = ch_versions.mix(NICHECOMPASS_ANALYSIS.out.versions.first())

    emit:
    nichecompass_dir = NICHECOMPASS_ANALYSIS.out.nichecompass_dir  // channel: [ val(meta), [ nichecompass_dir ] ]
    notebook         = NICHECOMPASS_ANALYSIS.out.notebook          // channel: [ val(meta), [ notebook.ipynb ] ]
    versions         = ch_versions                                 // channel: [ versions.yml ]
}
