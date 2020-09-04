/************************************************************************
* PREPARE CSS
*
* Extract centroids and cluster from cluster-files and prepare
* css files for newick utils.
************************************************************************/

process prepare_css {
  label 'preparecss'
  publishDir "${params.output}/${params.css_output}", mode: 'copy', pattern: '*map'

  input:
    path(clusterFile)

  output:
    tuple val("${clusterFile.baseName}"), path ("${clusterFile.baseName}_cluster_css.map"),  path ("${clusterFile.baseName}_ornaments.map"), emit: css_cluster
    //tuple val(clusterFile), path ("${clusterFile.baseName}_ornaments.map"), emit: css_ornaments

  script:
  """
  python3 ${projectDir}/bin/prepare_css.py "${clusterFile}"

  """

}