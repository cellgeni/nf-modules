include { SEGTRAQ_QC as SEGTRAQ_QC_MODULE } from '../../../modules/sanger-cellgeni/segtraq/qc/main'

workflow SEGTRAQ_QC {
    take:
    ch_xenium_bundles // channel: [ val(meta), path(xenium_bundle), val(metrics) ]

    main:
    SEGTRAQ_QC_MODULE(ch_xenium_bundles)
    ch_baseline_obs = SEGTRAQ_QC_MODULE.out.baseline_obs
    ch_baseline_summary = SEGTRAQ_QC_MODULE.out.baseline_summary
    ch_clustering_obs = SEGTRAQ_QC_MODULE.out.clustering_obs
    ch_clustering_metrics = SEGTRAQ_QC_MODULE.out.clustering_metrics
    ch_region_similarity_obs = SEGTRAQ_QC_MODULE.out.region_similarity_obs
    ch_volume_obs = SEGTRAQ_QC_MODULE.out.volume_obs
    ch_versions = SEGTRAQ_QC_MODULE.out.versions

    emit:
    baseline_obs          = ch_baseline_obs
    baseline_summary      = ch_baseline_summary
    clustering_obs        = ch_clustering_obs
    clustering_metrics    = ch_clustering_metrics
    region_similarity_obs = ch_region_similarity_obs
    volume_obs            = ch_volume_obs
    versions              = ch_versions
}
