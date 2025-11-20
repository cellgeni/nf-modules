process SPATIAL_NICHECOMPASSANALYSIS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/cellgeni/nichecompass:0.3.0"

    stageInMode 'copy'

    input:
    tuple val(meta), path(nichecompass_model), path(nichecompass_run_config), path(nichecompass_data)

    output:
    tuple val(meta), path("${out_dir}"), emit: nichecompass_dir
    tuple val(meta), path("${out_dir}_analysis.ipynb"), emit: notebook
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_nichecompass"
    """
    mkdir ${out_dir}
    papermill \\
        "${moduleDir}/resources/usr/bin/nichecompass_analyse_sample_integration.ipynb" \\
        "${out_dir}_analysis.ipynb" \\
        -p run_root "${out_dir}" \\
        -p cfg_path "${nichecompass_run_config}" \\
        -p model_folder_path "${nichecompass_model}" \\
        ${args} \\
        --kernel python3 \\
        --request-save-on-cell-execute \\
        --progress-bar \\
        --log-level INFO \\
        --log-output


  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
  END_VERSIONS
  """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_dir = "${prefix}_nichecompass"
    """
    mkdir -p "${out_dir}/figures"
    touch "${out_dir}_analysis.ipynb"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      nichecompass: \$(pip show nichecompass | grep Version | sed -e "s/Version: //g")
  END_VERSIONS
  """
}
