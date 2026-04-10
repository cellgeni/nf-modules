include { SEGGER_SEGMENT      } from '../../../modules/sanger-cellgeni/segger/segment/main'
include { SEGGER_EXPORT       } from '../../../modules/sanger-cellgeni/segger/export/main'
include { SEGGER_XENIUMBUNDLE } from '../../../modules/sanger-cellgeni/segger/xeniumbundle/main'

workflow SEGGER_PREDICT {
    take:
    ch_input // channel: [ val(meta), path(input_dir) ]

    main:
    def ch_versions = channel.empty()

    SEGGER_SEGMENT(ch_input)
    ch_versions = ch_versions.mix(SEGGER_SEGMENT.out.versions.first())

    // Join segmentation result with the original input_dir to use as source_dir for export
    def ch_export = SEGGER_SEGMENT.out.segmentation
        .join(ch_input, by: 0)

    SEGGER_EXPORT(ch_export)
    ch_versions = ch_versions.mix(SEGGER_EXPORT.out.versions.first())

    // Join export_dir with original input_dir (xenium bundle), reordering to [ meta, xenium_dir, export_dir ]
    def ch_bundle = SEGGER_EXPORT.out.export_dir
        .join(ch_input, by: 0)
        .map { meta, export_dir, xenium_dir -> [ meta, xenium_dir, export_dir ] }

    SEGGER_XENIUMBUNDLE(ch_bundle)
    ch_versions = ch_versions.mix(SEGGER_XENIUMBUNDLE.out.versions.first())

    emit:
    segmentation  = SEGGER_SEGMENT.out.segmentation        // channel: [ val(meta), path(parquet) ]
    export_dir    = SEGGER_EXPORT.out.export_dir           // channel: [ val(meta), path(export_dir) ]
    xenium_bundle = SEGGER_XENIUMBUNDLE.out.xenium_bundle // channel: [ val(meta), path(xenium_bundle) ]
    versions      = ch_versions                            // channel: [ versions.yml ]
}
