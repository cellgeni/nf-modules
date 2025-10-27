def metaToTsv(meta) {
    def tsv_string = meta
        .findAll { key, value -> key != 'id' && value }
        .collectMany { key, value ->
            value
                .toString()
                .split(/\s*,\s*/)
                .collect { it.trim() }
                .findAll { it }
                .collect { v -> "${key}\t${v}" }
        }
        .join('\\n')
        .stripIndent()
        .replaceAll('"', '\\\\"')
    // remove leading whitespace and escape quotes
    return tsv_string
}

process IRODS_ATTACHMETADATA {
    tag "Attaching metadata for ${meta.id}"

    input:
    tuple val(meta), val(irodspath)

    output:
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def irodspath = irodspath.replaceFirst('/$', '')
    def meta_tsv = metaToTsv(meta)
    def delimiter = task.ext.delimiter ?: ""
    """
    # Create tsv file with metadata
    set -euo pipefail
    echo -e "${meta_tsv}" > metadata.tsv

    # Check if irodspath exists
    if ils -d "${irodspath}" | grep -q ':\$'; then
        resource="-C"
    elif ils -d "${irodspath}"; then
        resource="-d"
    else
        echo "Error: iRODS path ${irodspath} does not exist."
        exit 1
    fi

    # Get existing metadata from iRODS
    get_metadata.sh \$resource "${irodspath}" > existing_metadata.csv

    echo "Existing metadata for ${irodspath}:"
    cat existing_metadata.csv

    # Remove existing metadata if specified
    if [ "${task.ext.remove_existing_metadata}" == "true" ]; then
        echo "Removing existing metadata for ${irodspath}"
        imeta rmw \$resource "${irodspath}" % % || echo "No metadata to remove (this is OK)"
        :> existing_metadata.csv # clear file
    fi

    # Load metadata to iRODS
    echo "Current metadata for ${irodspath}:"
    get_metadata.sh \$resource "${irodspath}"
    set +e
    while IFS=\$'\\t' read -r key value; do
        [[ -z "\$key" || -z "\$value" ]] && continue  # skip empty lines

        # Check if value contains semicolon delimiter
        if [[ -n "${delimiter}" && "\$value" == *"${delimiter}"* ]]; then
            # Split by semicolon and process each value separately
            IFS='${delimiter}' read -ra VALUES <<< "\$value"
            for val in "\${VALUES[@]}"; do
                val=\$(echo "\$val" | xargs)  # trim whitespace
                [[ -z "\$val" ]] && continue  # skip empty values
                
                # Check if the key value pair already exists in iRODS metadata
                if grep -qzP "\${key},\${val}" existing_metadata.csv; then
                    echo "[SKIP] \$key=\$val already present"
                else
                    echo "Adding \$key=\$val to iRODS metadata"
                    imeta add \$resource "${irodspath}" "\$key" "\$val"
                fi
            done
        else
            # Process single value as before
            if grep -qzP "\${key},\${value}" existing_metadata.csv; then
                echo "[SKIP] \$key=\$value already present"
            else
                echo "Adding \$key=\$value to iRODS metadata"
                imeta add \$resource "${irodspath}" "\$key" "\$value"
            fi
        fi
    done < metadata.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
    END_VERSIONS
    """
}
