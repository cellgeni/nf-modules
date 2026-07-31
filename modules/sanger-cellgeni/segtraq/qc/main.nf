process SEGTRAQ_QC {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/segtraq:0.0.4"

    input:
    tuple val(meta), path(xenium_dir), val(metrics)

    output:
    tuple val(meta), path("${prefix}_baseline_obs.csv"), emit: baseline_obs, optional: true
    tuple val(meta), path("${prefix}_baseline_summary.json"), emit: baseline_summary, optional: true
    tuple val(meta), path("${prefix}_clustering_stability_obs.csv"), emit: clustering_obs, optional: true
    tuple val(meta), path("${prefix}_clustering_stability_metrics.json"), emit: clustering_metrics, optional: true
    tuple val(meta), path("${prefix}_region_similarity_obs.csv"), emit: region_similarity_obs, optional: true
    tuple val(meta), path("${prefix}_volume_similarity_obs.csv"), emit: volume_obs, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${moduleDir}/resources/usr/bin/segtraq_qc.py \\
        --xenium_dir "${xenium_dir}" \\
        --prefix "${prefix}" \\
        --metrics "${metrics}" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segtraq: \$(uv pip show segtraq | grep "^Version" | sed 's/Version: //')
        segtraq_qc: \$(python3 ${moduleDir}/resources/usr/bin/segtraq_qc.py --version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def metricList = metrics instanceof List
        ? metrics.collect { metric -> "${metric}".trim().toLowerCase() }
        : "${metrics}".split(',').collect { metric -> metric.trim().toLowerCase() }
    if (metricList.contains('all')) {
        metricList = ['baseline', 'clustering_stability', 'region_similarity', 'volume']
    }
    """
    if ${metricList.contains('baseline')}; then
        touch "${prefix}_baseline_obs.csv"
        echo '{"num_cells":0,"perc_unassigned_transcripts":0.0}' > "${prefix}_baseline_summary.json"
    fi
    if ${metricList.contains('clustering_stability') || metricList.contains('clusteringstability') || metricList.contains('clustering')}; then
        touch "${prefix}_clustering_stability_obs.csv"
        echo '{}' > "${prefix}_clustering_stability_metrics.json"
    fi
    if ${metricList.contains('region_similarity') || metricList.contains('regionsimilarity') || metricList.contains('region')}; then
        touch "${prefix}_region_similarity_obs.csv"
    fi
    if ${metricList.contains('volume')}; then
        touch "${prefix}_volume_similarity_obs.csv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segtraq: "stub"
        segtraq_qc: "stub"
    END_VERSIONS
    """
}
