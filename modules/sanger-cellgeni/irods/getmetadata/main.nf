process IRODS_GETMETADATA {
    tag "Getting metadata for ${irodspath}"

    input:
    tuple val(meta), val(irodspath)

    output:
    tuple val(meta), path("irods_metadata.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    # Check if irodspath exists
    if ils -d "${irodspath}" | grep -q ':\$'; then
        resource="-C"
    elif ils -d "${irodspath}"; then
        resource="-d"
    else
        echo "Error: iRODS path ${irodspath} does not exist."
        exit 1
    fi

    # Get metadata from iRODS
    imeta ls \$resource ${irodspath} \
        | (grep -E 'attribute|value|units' || true) \
        | sed -e 's/^attribute: //' -e 's/^value: //' -e 's/^units: //' \
        | sed -e "s/\\\"/'/g" \
        | awk 'NR%3!=0 {printf "\\\"%s\\\",", \$0} NR%3==0 {printf "\\\"%s\\\"\\n", \$0}' > irods_metadata.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
        awk: \$(awk --version | head -n1)
        sed: \$(sed --version | head -n1 | awk '{ print \$4 }')
    END_VERSIONS
    """

    stub:
    """
    touch irods_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
        awk: \$(awk --version | head -n1)
        sed: \$(sed --version | head -n1 | awk '{ print \$4 }')
    END_VERSIONS
    """
}
