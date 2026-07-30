include { SEGTRAQ_QC as SEGTRAQ_QC_MODULE } from '../../../modules/sanger-cellgeni/segtraq/qc/main'

workflow SEGTRAQ_QC {
    take:
    ch_xenium_bundles // channel: [ val(meta), path(xenium_bundle), val(metrics) ]

    main:
    SEGTRAQ_QC_MODULE(ch_xenium_bundles)

    emit:
    baseline_obs          = SEGTRAQ_QC_MODULE.out.baseline_obs
    baseline_summary      = SEGTRAQ_QC_MODULE.out.baseline_summary
    clustering_obs        = SEGTRAQ_QC_MODULE.out.clustering_obs
    clustering_metrics    = SEGTRAQ_QC_MODULE.out.clustering_metrics
    region_similarity_obs = SEGTRAQ_QC_MODULE.out.region_similarity_obs
    volume_obs            = SEGTRAQ_QC_MODULE.out.volume_obs
    versions              = SEGTRAQ_QC_MODULE.out.versions
}
