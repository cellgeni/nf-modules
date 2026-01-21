include { SPATIAL_TERRA } from '../../../modules/sanger-cellgeni/spatial/terra/main'
include { SPATIAL_CONCATH5ADS } from '../../../modules/sanger-cellgeni/spatial/concath5ads/main'

workflow TERRA {
    take:
    ch_anndatas // channel: [ val(meta), [ anndata ] ]

    main:

    ch_versions = channel.empty()

    SPATIAL_TERRA(ch_anndatas)
    ch_versions = ch_versions.mix(SPATIAL_TERRA.out.versions.first())

    SPATIAL_CONCATH5ADS(SPATIAL_TERRA.out.processed_anndata.collect())
    ch_versions = ch_versions.mix(SPATIAL_CONCATH5ADS.out.versions.first())

    emit:
    processed_anndata = SPATIAL_TERRA.out.processed_anndata // channel: [ val(meta), [ processed_anndata ] ]
    concatenated_h5ad = SPATIAL_CONCATH5ADS.out.concatenated_h5ad // channel: [ val(meta), [ concatenated_h5ad ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
