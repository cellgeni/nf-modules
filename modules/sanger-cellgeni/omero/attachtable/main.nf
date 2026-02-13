process OMERO_ATTACHTABLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE'
        : 'biocontainers/YOUR-TOOL-HERE'}"

    input:
    tuple val(meta), path(csv), val(image_id), val(roi_id), val(host), val(port)

    secret 'OMERO_USER'
    secret 'OMERO_PASSWORD'

    output:
    tuple val(meta), path("*.ann.csv"), emit: ann_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    attach_csv.py \\
        --host ${host} \\
        --port ${port} \\
        --user \$OMERO_USER \\
        --password \$OMERO_PASSWORD \\
        --image-id ${image_id} \\
        --roi-id ${roi_id} \\
        --csv ${csv} \\
        --ann-csv ${prefix}.ann.csv \\
        ${args}

    python - <<'PY' > versions.yml
    import importlib
    import platform

    def get_version(name):
        try:
            mod = importlib.import_module(name)
            return getattr(mod, "__version__", "unknown")
        except Exception:
            return "unknown"

    print("${task.process}:")
    print(f"  python: {platform.python_version()}")
    print(f"  omero2pandas: {get_version('omero2pandas')}")
    PY
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    cat <<-EOF > ${prefix}.ann.csv
    ann_id,image_id,roi_id
    0,${image_id},${roi_id}
    EOF

        cat <<-EOF > versions.yml
    ${task.process}:
    python: stub
    omero2pandas: stub
    EOF
    """
}
