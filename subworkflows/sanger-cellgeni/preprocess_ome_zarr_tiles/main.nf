include { BIOFORMATS2RAWCOMPANION      } from "../../../modules/sanger-cellgeni/bioformats2rawcompanion"
include { IMAGING_GENERATECOMPANION    } from "../../../modules/sanger-cellgeni/imaging/generatecompanion"
include { OMEZARR_CREATEEMPTYCONTAINER } from "../../../modules/sanger-cellgeni/omezarr/createemptycontainer"
include { OMEZARR_PREPROCESS           } from "../../../modules/sanger-cellgeni/omezarr/preprocess"
include { OMEZARR_VALIDATE             } from "../../../modules/sanger-cellgeni/omezarr/validate"

workflow PREPROCESS_OME_ZARR_TILES {
    take:
    master_file_ch // channel: [ val(meta), val(round), path(image_root_folder), val(master_file_name) ]
    psf_folder // channel: [ path(psf_folder) ]

    main:
    ch_versions = channel.empty()

    // 1. Convert the raw tiled acquisition into a single OME-Zarr container/plate,
    //    unless image_root_folder is already an OME-Zarr store (a directory whose
    //    name ends in ".zarr"), in which case the conversion is skipped entirely
    //    and that folder is used directly for the steps below.
    master_file_branches = master_file_ch.branch { meta, _round, image_root_folder, _master_file_name ->
        already_zarr: image_root_folder.isDirectory() && image_root_folder.name.endsWith('.zarr')
        raw: true
    }

    BIOFORMATS2RAWCOMPANION(
        master_file_branches.raw.map { meta, _round, image_root_folder, master_file_name ->
            [meta, image_root_folder.resolve(master_file_name), image_root_folder]
        }
    )
    ch_versions = ch_versions.mix(BIOFORMATS2RAWCOMPANION.out.versions_bioformats2rawcompanion.first())

    ch_already_zarr = master_file_branches.already_zarr.map { meta, _round, image_root_folder, _master_file_name -> [meta, image_root_folder] }

    ch_input_ome_zarr = BIOFORMATS2RAWCOMPANION.out.ome_zarr.mix(ch_already_zarr)
    // [ meta, ome_zarr ]

    // Tile manifest: which (image_id, hcs_path) tiles exist for this dataset.
    IMAGING_GENERATECOMPANION(master_file_ch)
    ch_versions = ch_versions.mix(IMAGING_GENERATECOMPANION.out.versions.first())

    tiles = IMAGING_GENERATECOMPANION.out.csv
        .splitCsv(header: true, strip: true, sep: ",")
        .map { meta, row -> [meta, row.image_id, row.path] }
    // [ meta, image_id, hcs_path ]

    tiles_per_dataset = tiles
        .map { meta, image_id, hcs_path -> [meta, [image_id, hcs_path]] }
        .groupTuple(by: 0)
        .map { meta, tile_list -> [meta, groovy.json.JsonOutput.toJson(tile_list)] }
    // [ meta, tiles_json ]

    // 2. Create the output OME-Zarr store upfront, with one correctly-shaped empty
    //    image slot pre-registered per tile, so step 3 can write into it safely in
    //    parallel without racing on shared zarr group/plate metadata.
    OMEZARR_CREATEEMPTYCONTAINER(
        tiles_per_dataset.combine(ch_input_ome_zarr, by: 0)
    )
    ch_versions = ch_versions.mix(OMEZARR_CREATEEMPTYCONTAINER.out.versions_omezarr_createemptycontainer.first())
    ch_out_zarr = OMEZARR_CREATEEMPTYCONTAINER.out.out_zarr
    // [ meta, out_zarr ]

    // 3. Process every tile in parallel, writing straight into its pre-created slot.
    tiles_to_process = tiles
        .combine(ch_input_ome_zarr, by: 0)
        .combine(ch_out_zarr, by: 0)
        .map { meta, image_id, hcs_path, root_folder, out_zarr ->
            [meta, root_folder, image_id, hcs_path, out_zarr]
        }
    OMEZARR_PREPROCESS(tiles_to_process, psf_folder)
    ch_versions = ch_versions.mix(OMEZARR_PREPROCESS.out.versions_omezarr_preprocess.first())

    // 4. Once every tile for a dataset has been written, validate the finished store.
    //    groupTuple only flushes a key once the whole upstream channel completes, so
    //    all OMEZARR_PREPROCESS tasks across the run finish before any validation starts.
    completed_tiles = OMEZARR_PREPROCESS.out.fovs
        .map { meta, out_zarr, image_id, hcs_path -> [meta, out_zarr, [image_id, hcs_path]] }
        .groupTuple(by: [0, 1])
        .map { meta, out_zarr, tile_list -> [meta, out_zarr, groovy.json.JsonOutput.toJson(tile_list)] }
    OMEZARR_VALIDATE(completed_tiles)
    ch_versions = ch_versions.mix(OMEZARR_VALIDATE.out.versions_omezarr_validate.first())

    emit:
    companion_tiles = IMAGING_GENERATECOMPANION.out.companion.join(OMEZARR_VALIDATE.out.out_zarr, by: 0) // [ meta, companion, out_zarr ]
    versions        = ch_versions // channel: [ versions.yml ]
}
