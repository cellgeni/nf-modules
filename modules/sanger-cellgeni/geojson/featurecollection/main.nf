process GEOJSON_FEATURECOLLECTION {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(geojson)

    output:
    tuple val(meta), path("${prefix}.geojson"), emit: geojson
    eval("python3 ${moduleDir}/resources/usr/bin/convert_geojson_featurecollection.py --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    """
    python3 ${moduleDir}/resources/usr/bin/convert_geojson_featurecollection.py \
        --input $geojson \
        --output ${prefix}.geojson \
        $args
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    cat << 'JSON' > ${prefix}.geojson
    {
      "type": "FeatureCollection",
      "features": []
    }
    JSON
    """
}
