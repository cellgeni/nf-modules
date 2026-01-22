process CSV_CONCAT {
    tag "Concatenating CSV files"
    container 'docker://quay.io/cellgeni/metacells-python:latest'

    input:
    tuple val(meta), path(csv_files, name: "input/*.csv")

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "concatenated"
    def args = task.ext.args ?: ""
    """
    concat.py --input ${csv_files} --prefix "${prefix}" ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | awk '{print \$2}')
        pandas: \$( python -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """
}
