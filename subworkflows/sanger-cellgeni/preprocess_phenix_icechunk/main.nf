include { PHENIXICECHUNK         } from "../../../modules/sanger-cellgeni/phenixicechunk"
include { MULTIVIEWSTITCHER      } from "../../../modules/sanger-cellgeni/multiviewstitcher"
include { ICECHUNK_STITCHEDPLATE } from "../../../modules/sanger-cellgeni/icechunk/stitchedplate"

workflow PREPROCESS_PHENIX_ICECHUNK {
    take:
    measurements_ch // channel: [ val(meta), path(image_dir), val(index_name) ], meta: [ id, round_index ]

    main:
    // One repo per measurement: make meta.id unique per round and keep the original id
    // as meta.sample to group the rounds back together
    measurements = measurements_ch.map { meta, image_dir, index_name ->
        [meta + [id: "${meta.id}_r${meta.round_index}".toString(), sample: meta.id], image_dir, index_name]
    }
    PHENIXICECHUNK(measurements)

    // One task per round and well, so the wells are stitched in parallel. Rounds are not
    // registered to each other. stitch_wells.py gets the well and keeps the 3D stacks
    // through ext.args: { "--projection none --wells ${meta.well}" }
    well_repos = PHENIXICECHUNK.out.repo
        .join(PHENIXICECHUNK.out.wells)
        .flatMap { meta, repo, wells_txt ->
            wells_txt.readLines().collect { line -> line.trim() }.findAll { well -> well }.collect { well ->
                [meta + [id: "${meta.id}_${well}".toString(), well: well], repo]
            }
        }
    MULTIVIEWSTITCHER(well_repos)

    // Collect every stitched well of a sample, with the source repos of all its rounds,
    // into one plate, each round an acquisition. Sorted so the task inputs, and therefore
    // -resume, are stable.
    source_repos = PHENIXICECHUNK.out.repo
        .map { meta, repo -> [meta.sample, meta.round_index, repo] }
        .groupTuple(by: 0)
        .map { sample, rounds, repos ->
            def ordered = [rounds, repos].transpose().sort { a, b -> (a[0] as Integer) <=> (b[0] as Integer) }
            [sample, ordered.collect { round_repo -> round_repo[1] }, ordered.collect { round_repo -> round_repo[0] }]
        }
    plate_input = MULTIVIEWSTITCHER.out.ome_zarr
        .map { meta, ome_zarr -> [meta.sample, ome_zarr, [well: meta.well, rounds: [meta.round_index]]] }
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
    stitched = MULTIVIEWSTITCHER.out.ome_zarr // channel: [ val(meta), path(ome.zarr) ], 3D stack per round and well
    qc       = MULTIVIEWSTITCHER.out.qc // channel: [ val(meta), path(qc.png) ], per round and well
    plate    = ICECHUNK_STITCHEDPLATE.out.repo // channel: [ val(meta), path(stitched.icechunk) ], one HCS plate per sample
}
