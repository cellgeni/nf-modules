include { PHENIXICECHUNK         } from "../../../modules/sanger-cellgeni/phenixicechunk"
include { ASHLAR                 } from "../../../modules/nf-core/ashlar"
include { MULTIVIEWSTITCHER      } from "../../../modules/sanger-cellgeni/multiviewstitcher"
include { ICECHUNK_STITCHEDPLATE } from "../../../modules/sanger-cellgeni/icechunk/stitchedplate"

workflow PREPROCESS_PHENIX_ICECHUNK {
    take:
    measurements_ch // channel: [ val(meta), path(image_dir), val(index_name) ], meta: [ id, round_index ]
    stitch_3d // boolean: true = stitch the 3D stacks of every round and well with multiview-stitcher, false = stitch and register every well across rounds with ASHLAR
    dfp_folder // ASHLAR dark field profile(s), or []
    ffp_folder // ASHLAR flat field profile(s), or []

    main:
    // One repo per measurement: make meta.id unique per round and keep the original id
    // as meta.sample to group the rounds back together
    measurements = measurements_ch.map { meta, image_dir, index_name ->
        [meta + [id: "${meta.id}_r${meta.round_index}".toString(), sample: meta.id], image_dir, index_name]
    }
    PHENIXICECHUNK(measurements)

    // 2D: split per well and order each well's measurements by round_index, so ashlar
    // registers every later cycle against the first one. Ashlar needs z projected fields,
    // so PHENIXICECHUNK must export <plate>_<well>_z<method>.ome.tiff (ext.args `--zproject max --ome-tiff .`)
    well_cycles = PHENIXICECHUNK.out.ome_tiff
        .flatMap { meta, tiffs ->
            (tiffs instanceof List ? tiffs : [tiffs]).collect { tiff -> [meta, tiff] }
        }
        .map { meta, tiff ->
            def well = (tiff.name =~ /_([A-Z]+\d+)(?:_z[a-z]+)?\.ome\.tiff$/)[0][1]
            [[id: "${meta.sample}_${well}".toString(), sample: meta.sample, well: well], meta.round_index, tiff]
        }
        .groupTuple(by: 0)
        .map { meta, rounds, tiffs ->
            def ordered = [rounds, tiffs].transpose().sort { a, b -> (a[0] as Integer) <=> (b[0] as Integer) }
            [meta + [rounds: ordered.collect { round_tiff -> round_tiff[0] }], ordered.collect { round_tiff -> round_tiff[1] }]
        }

    // 3D: one task per round and well, so the wells are stitched in parallel. stitch_wells.py
    // gets the well and keeps the stacks through ext.args: { "--projection none --wells ${meta.well}" }
    well_repos = PHENIXICECHUNK.out.repo
        .join(PHENIXICECHUNK.out.wells)
        .flatMap { meta, repo, wells_txt ->
            wells_txt.readLines().collect { line -> line.trim() }.findAll { well -> well }.collect { well ->
                [meta + [id: "${meta.id}_${well}".toString(), well: well], repo]
            }
        }

    if (stitch_3d) {
        MULTIVIEWSTITCHER(well_repos)
    }
    else {
        ASHLAR(well_cycles, dfp_folder, ffp_folder)
    }
    // [ meta, image, [ well, rounds ] ]: rounds are the round_index values whose channels the image holds
    stitched = stitch_3d
        ? MULTIVIEWSTITCHER.out.ome_zarr.map { meta, ome_zarr -> [meta, ome_zarr, [well: meta.well, rounds: [meta.round_index]]] }
        : ASHLAR.out.tif.map { meta, tif -> [meta, tif, [well: meta.well, rounds: meta.rounds]] }

    // Collect every stitched well of a sample, with the source repos of all its rounds,
    // into one plate. Sorted so the task inputs, and therefore -resume, are stable.
    source_repos = PHENIXICECHUNK.out.repo
        .map { meta, repo -> [meta.sample, meta.round_index, repo] }
        .groupTuple(by: 0)
        .map { sample, rounds, repos ->
            def ordered = [rounds, repos].transpose().sort { a, b -> (a[0] as Integer) <=> (b[0] as Integer) }
            [sample, ordered.collect { round_repo -> round_repo[1] }, ordered.collect { round_repo -> round_repo[0] }]
        }
    plate_input = stitched
        .map { meta, image, entry -> [meta.sample, image, entry] }
        .groupTuple(by: 0)
        .map { sample, images, layout ->
            def ordered = [images, layout].transpose().sort { a, b -> "${a[1].well} ${a[1].rounds}".toString() <=> "${b[1].well} ${b[1].rounds}".toString() }
            [sample, ordered.collect { image_entry -> image_entry[0] }, ordered.collect { image_entry -> image_entry[1] }]
        }
        .join(source_repos)
        .map { sample, images, layout, repos, rounds -> [[id: sample], images, layout, repos, rounds] }
    ICECHUNK_STITCHEDPLATE(plate_input)

    emit:
    repo     = PHENIXICECHUNK.out.repo // channel: [ val(meta), path(repo.icechunk) ], one per measurement
    cycles   = stitch_3d ? channel.empty() : well_cycles // channel: [ val(meta), [ path(ome.tiff) ] ], one per well ordered by round_index; empty when stitch_3d
    stitched = stitched.map { meta, image, _entry -> [meta, image] } // channel: [ val(meta), path(ome.tif | ome.zarr) ], one per well (2D) or per round and well (3D)
    qc       = stitch_3d ? MULTIVIEWSTITCHER.out.qc : channel.empty() // channel: [ val(meta), path(qc.png) ], 3D only
    plate    = ICECHUNK_STITCHEDPLATE.out.repo // channel: [ val(meta), path(stitched.icechunk) ], one HCS plate per sample
}
