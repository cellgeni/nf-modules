#!/usr/bin/env nextflow

include { SPATIAL_SCRAFTPALOM } from '../../../modules/sanger-cellgeni/spatial/scraftpalom/main'
include { OMETIF_CONCATENATE  } from '../../../modules/sanger-cellgeni/ometif/concatenate/main'

workflow ALIGN_SERIES_TO_REFERENCE {
    take:
    ch_images // channel: [ val(meta), path(images) ] ordered list of images to align; one image is used unchanged as the reference

    main:
    ch_versions = channel.empty()

    method     = params.alignment_method ?: 'scraftpalom'
    ref_index  = params.align_reference_index ?: 0

    // Split the ordered image list into a single reference image and the
    // remaining moving images, keyed by ref_index (defaults to the first image).
    ch_ref_moving = ch_images.map { meta, images ->
        def imgs = images instanceof List ? images : [images]
        def ref = imgs[ref_index]
        def moving = imgs.indices.findAll { it != ref_index }.collect { imgs[it] }
        tuple(meta, ref, moving)
    }

    // One registration task per moving image, all sharing the same reference.
    // Each task's moving-list position is stashed in the meta map (removed again
    // right after) so the aligned images can be put back in their original order
    // afterwards -- registration tasks can finish in any order.
    ch_pairs = ch_ref_moving.flatMap { meta, ref, moving ->
        moving.withIndex().collect { m, idx -> tuple(meta + [_align_order: idx], ref, m) }
    }

    if (method == 'scraftpalom') {
        SPATIAL_SCRAFTPALOM(ch_pairs)
        ch_aligned_tagged = SPATIAL_SCRAFTPALOM.out.aligned_moving_image // [ tagged_meta, aligned_image ]
        ch_versions = ch_versions.mix(SPATIAL_SCRAFTPALOM.out.versions.first())
    } else {
        error "Unsupported alignment method '${method}'. Currently supported methods: ['scraftpalom']"
    }

    ch_aligned_ordered = ch_aligned_tagged.map { tagged_meta, img ->
        def order = tagged_meta._align_order
        def meta = tagged_meta.findAll { k, _v -> k != '_align_order' }
        tuple(meta, order, img)
    }

    // Reassemble reference + aligned images into a single ordered series per
    // sample: the reference first, then the aligned images in their original
    // (pre-alignment) relative order.
    ch_aligned_series = ch_aligned_ordered
        .map { meta, order, img -> tuple(meta, [order, img]) }
        .groupTuple(by: 0)
        .map { meta, order_img_pairs -> tuple(meta, order_img_pairs.sort { it[0] }.collect { it[1] }) }
        .join(ch_ref_moving.map { meta, ref, _moving -> tuple(meta, ref) }, by: 0)
        .map { meta, aligned_images, ref -> tuple(meta, [ref] + aligned_images) }

    OMETIF_CONCATENATE(ch_aligned_series)
    ch_versions = ch_versions.mix(OMETIF_CONCATENATE.out.versions.first())

    emit:
    reference  = ch_ref_moving.map { meta, ref, _moving -> tuple(meta, ref) } // channel: [ val(meta), path(reference_image) ]
    aligned    = ch_aligned_ordered.map { meta, _order, img -> tuple(meta, img) } // channel: [ val(meta), path(aligned_image) ] one entry per moving image
    hyperstack = OMETIF_CONCATENATE.out.hyperstack // channel: [ val(meta), path(hyperstack.ome.tif) ] reference + aligned images concatenated in original series order, channel names prefixed with their 1-based cycle position
    versions   = ch_versions // channel: [ versions.yml ]
}
