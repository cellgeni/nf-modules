include { SPATIALDATA_EXPORT             } from '../../../modules/sanger-cellgeni/spatialdata/export/main'
include { SPATIAL_GENERATEVITESSCECONFIG } from '../../../modules/sanger-cellgeni/spatial/generatevitessceconfig/main'

workflow SPATIALDATA2VITESSCE {
    take:
    ch_input // channel: [ val(meta), path(bundle) ]
    raw_name // val: image name under images/
    label_name // val: label name under labels/
    table_name // val: table name under tables/
    http_base_url // val: base HTTP URL for vitessce to serve the zarr

    main:
    ch_versions = channel.empty()

    SPATIALDATA_EXPORT(ch_input)
    ch_versions = ch_versions.mix(SPATIALDATA_EXPORT.out.versions.first())

    SPATIAL_GENERATEVITESSCECONFIG(
        SPATIALDATA_EXPORT.out.zarr.combine(raw_name).combine(label_name).combine(table_name).combine(http_base_url)
    )
    ch_versions = ch_versions.mix(SPATIAL_GENERATEVITESSCECONFIG.out.versions.first())

    emit:
    zarr            = SPATIALDATA_EXPORT.out.zarr // channel: [ val(meta), path(*.zarr) ]
    vitessce_config = SPATIAL_GENERATEVITESSCECONFIG.out.vitessce_config // channel: [ val(meta), path(*_vitessce_config.json) ]
    versions        = ch_versions // channel: [ versions.yml ]
}
