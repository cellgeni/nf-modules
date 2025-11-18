process IRODS_STOREFILE {
    tag "Loading ${irodspath}"

    input:
    tuple val(meta), path(file), val(irodspath)

    output:
    tuple val(meta), val(irodspath), env('md5'), env('irods_md5'), emit: md5
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: "-KV -f -X restart.txt --retries 10 --acl 'read public#archive'"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    module load cellgen/irods

    # calculate MD5
    md5=\$(md5sum "${file}" | awk '{print \$1}')

    # create iRODS directory if it doesn't exist
    irodsdir=\$(dirname "${irodspath}")
    imkdir -p "\$irodsdir"

    # Load file to iRODS
    echo "Loading ${file} to iRODS at ${irodspath}"
    iput ${args} \
        -N ${task.cpus} \
        --metadata="md5;\${md5};;" \
        "${file}" "${irodspath}"

    # Calculate iRODS md5
    sleep 1 # wait for iRODS to do it's thing
    irods_md5=\$(ichksum "${irodspath}" | awk '{print \$NF}')

    # Compare iRODS md5 with local md5
    if [ "\$md5" != "\$irods_md5" ]; then
        echo "MD5 mismatch for ${file}: local \$md5, iRODS \$irods_md5"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
        awk: \$(awk --version | head -n1)
        md5sum: \$(md5sum --version | head -n1 | awk '{ print \$4 }')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: "-K -f -X restart.txt --retries 10"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # calculate MD5
    md5=\$(md5sum "${file}" | awk '{print \$1}')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        irods: \$(ienv | grep version | awk '{ print \$3 }')
        awk: \$(awk --version | head -n1)
        md5sum: \$(md5sum --version | head -n1 | awk '{ print \$4 }')
    END_VERSIONS
    """
}
