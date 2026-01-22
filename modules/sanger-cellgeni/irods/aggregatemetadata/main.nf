process IRODS_AGGREGATEMETADATA {
    tag "Aggregating metadata for ${meta.id}"
    container 'docker://quay.io/cellgeni/metacells-python:latest'

    input:
    tuple val(meta), path(irods_metadata, name: "input.csv")

    output:
    tuple val(meta), path("metadata.csv"), emit: csv
    tuple val(meta), path("metadata.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: '--dup-sep ";" --index_name "id"'
    """
    aggregate_metadata.py \
        ${args} \
        --input input.csv \
        --id ${meta.id}
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '--dup-sep "," --index_name "id"'
    """
    touch metadata.csv metadata.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """
}
