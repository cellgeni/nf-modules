include { SEGTRAQ_BASELINE           } from '../../../modules/sanger-cellgeni/segtraq/baseline/main'
include { SEGTRAQ_CLUSTERINGSTABILITY } from '../../../modules/sanger-cellgeni/segtraq/clusteringstability/main'
include { SEGTRAQ_REGIONSIMILARITY    } from '../../../modules/sanger-cellgeni/segtraq/regionsimilarity/main'
include { SEGTRAQ_VOLUME              } from '../../../modules/sanger-cellgeni/segtraq/volume/main'

workflow SEGTRAQ_QC {
    take:
    ch_zarr // channel: [ val(meta), path(zarr_dir) ]

    main:
    def ch_versions = channel.empty()

    // All four modules are independent and run in parallel on the same zarr_dir
    SEGTRAQ_BASELINE(ch_zarr)
    ch_versions = ch_versions.mix(SEGTRAQ_BASELINE.out.versions.first())

    SEGTRAQ_CLUSTERINGSTABILITY(ch_zarr)
    ch_versions = ch_versions.mix(SEGTRAQ_CLUSTERINGSTABILITY.out.versions.first())

    SEGTRAQ_REGIONSIMILARITY(ch_zarr)
    ch_versions = ch_versions.mix(SEGTRAQ_REGIONSIMILARITY.out.versions.first())

    SEGTRAQ_VOLUME(ch_zarr)
    ch_versions = ch_versions.mix(SEGTRAQ_VOLUME.out.versions.first())

    emit:
    baseline_obs          = SEGTRAQ_BASELINE.out.obs           // channel: [ val(meta), path(*_baseline_obs.csv) ]
    baseline_summary      = SEGTRAQ_BASELINE.out.summary       // channel: [ val(meta), path(*_baseline_summary.json) ]
    clustering_obs        = SEGTRAQ_CLUSTERINGSTABILITY.out.obs     // channel: [ val(meta), path(*_clustering_stability_obs.csv) ]
    clustering_metrics    = SEGTRAQ_CLUSTERINGSTABILITY.out.metrics // channel: [ val(meta), path(*_clustering_stability_metrics.json) ]
    region_similarity_obs = SEGTRAQ_REGIONSIMILARITY.out.obs   // channel: [ val(meta), path(*_region_similarity_obs.csv) ]
    volume_obs            = SEGTRAQ_VOLUME.out.obs             // channel: [ val(meta), path(*_volume_similarity_obs.csv) ]
    versions              = ch_versions                        // channel: [ versions.yml ]
}
