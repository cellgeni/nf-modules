include { NICHECOMPASS_PREPROCESS      } from '../../../modules/sanger-cellgeni/nichecompass/preprocess/main'
include { NICHECOMPASS_TRAIN           } from '../../../modules/sanger-cellgeni/nichecompass/train/main'
include { SPATIAL_NICHECOMPASSANALYSIS } from '../../../modules/sanger-cellgeni/spatial/nichecompassanalysis/main'
include { SCRAFT_SPLITH5AD             } from '../../../modules/sanger-cellgeni/scraft/splith5ad/main'

workflow NICHECOMPASS {
    take:
    ch_input // channel: [ val(meta), path(h5ad), val(cell_type_key), val(sample_key), val(counts_key), val(spatial_key) ]

    main:
    def ch_versions = channel.empty()

    // Carry the four mandatory AnnData keys alongside the split step;
    // sample_key is reused as the obs column to split on
    def ch_keys = ch_input.map { meta, h5ad, cell_type_key, sample_key, counts_key, spatial_key ->
        [meta, cell_type_key, sample_key, counts_key, spatial_key]
    }

    SCRAFT_SPLITH5AD(
        ch_input.map { meta, h5ad, cell_type_key, sample_key, counts_key, spatial_key ->
            [meta, h5ad, sample_key]
        }
    )
    ch_versions = ch_versions.mix(SCRAFT_SPLITH5AD.out.versions.first())

    NICHECOMPASS_PREPROCESS(
        SCRAFT_SPLITH5AD.out.h5ad_splits.join(ch_keys)
    )
    ch_versions = ch_versions.mix(NICHECOMPASS_PREPROCESS.out.versions.first())

    NICHECOMPASS_TRAIN(NICHECOMPASS_PREPROCESS.out.preprocessed)
    ch_versions = ch_versions.mix(NICHECOMPASS_TRAIN.out.versions.first())

    SPATIAL_NICHECOMPASSANALYSIS(NICHECOMPASS_TRAIN.out.nichecompass_model)
    ch_versions = ch_versions.mix(SPATIAL_NICHECOMPASSANALYSIS.out.versions.first())

    emit:
    nichecompass_dir = SPATIAL_NICHECOMPASSANALYSIS.out.nichecompass_dir // channel: [ val(meta), path(nichecompass_dir) ]
    notebook         = SPATIAL_NICHECOMPASSANALYSIS.out.notebook // channel: [ val(meta), path(notebook.ipynb) ]
    versions         = ch_versions // channel: [ versions.yml ]
}
