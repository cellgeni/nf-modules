include { PE2OMETIF } from "../../../modules/sanger-cellgeni/pe2ometif/main"
include { IMAGING_ASHLARCOMPANION } from "../../../modules/sanger-cellgeni/imaging/ashlarcompanion/main"
include { IMAGING_GENERATECOMPANIONFROMFILES } from "../../../modules/sanger-cellgeni/imaging/generatecompanionfromfiles/main"

workflow PREPROCESS_TIFF_TILES_ASHLAR {
    take:
    ch_input    // channel: [ val(meta{id, round_index}), path(image_dir), val(index_name) ]
    dfp_folder  // path or []
    ffp_folder  // path or []
    is_plate    // boolean

    main:
    ch_versions = Channel.empty()

    // 1. Convert PerkinElmer tiles to OME-TIFFs and per-well companion files
    PE2OMETIF(ch_input)

    // 2. Flatten per-well companions; group by (id, well) across rounds in order
    per_well_grouped = PE2OMETIF.out.per_well_companions
        .transpose()
        .map { meta, companion ->
            def well = (companion.name =~ /_([A-Z]{1,2}\d+)\.companion\.ome$/)[0][1]
            [ [id: meta.id, well: well], meta.round_index, companion ]
        }
        .groupTuple(by: 0)
        .map { meta, rounds, companions ->
            def sorted = [rounds, companions].transpose().sort { a, b -> a[0] <=> b[0] }
            [ meta, sorted.collect { it[1] } ]
        }

    // 3. Collect all tif files per sample across rounds, in round order
    tifs_by_sample = PE2OMETIF.out.ome_tif
        .map { meta, tifs ->
            [ [id: meta.id], meta.round_index, tifs instanceof List ? tifs : [tifs] ]
        }
        .groupTuple(by: 0)
        .map { meta, rounds, tifs_per_round ->
            def sorted = [rounds, tifs_per_round].transpose().sort { a, b -> a[0] <=> b[0] }
            [ meta, sorted.collect { it[1] }.flatten() ]
        }

    // 4. Pair each well's ordered companions with all tifs for that sample
    //    ASHLAR reads the per-well companion to determine which tiles belong to each well.
    multi_cycle_well = per_well_grouped
        .map { meta, companions -> [ [id: meta.id], meta.well, companions ] }
        .combine(tifs_by_sample, by: 0)
        .map { id_meta, well, companions, tifs ->
            [
                [ id: "${id_meta.id}_${well}", sample_id: id_meta.id, well: well ],
                companions,
                tifs
            ]
        }

    // 5. Stitch and register each well across all cycles with ASHLAR
    IMAGING_ASHLARCOMPANION(multi_cycle_well, dfp_folder, ffp_folder, is_plate)
    ch_versions = ch_versions.mix(IMAGING_ASHLARCOMPANION.out.versions.first())

    // 6. Collect all per-well stitched images per sample; generate master companion
    stitched_by_sample = IMAGING_ASHLARCOMPANION.out.tif
        .map { meta, tif -> [ [id: meta.sample_id], tif ] }
        .groupTuple()

    IMAGING_GENERATECOMPANIONFROMFILES(
        stitched_by_sample.combine(channel.of("*.ome.tif"))
    )
    ch_versions = ch_versions.mix(IMAGING_GENERATECOMPANIONFROMFILES.out.versions.first())

    emit:
    stitched_tif = IMAGING_ASHLARCOMPANION.out.tif                  // channel: [ val(meta{id, sample_id, well}), path(tif) ]
    companion    = IMAGING_GENERATECOMPANIONFROMFILES.out.companion  // channel: [ val(meta{id=sample_id}), path(companion) ]
    versions     = ch_versions                                       // channel: [ path(versions.yml) ]
}
