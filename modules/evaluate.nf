/************************************************************************
* EVALUATE
*
* Calculate basic statistics for each cluster algorithm
************************************************************************/

process evaluate_cluster {
  label 'evaluate'
  publishDir "${params.output}/${params.eval_output}", mode: 'copy', pattern: "*txt"

  input:
    tuple val(name), path(clusterFile), path(newick), path(sequences), val(flag), val(addParams)

  output:
    path "${name}_stats.out", emit: eval_result
    path ("${name}_taxonomy_info.txt") optional true

  script:
  def flag = flag == false ? '' : flag
  def addParams = addParams == false ? '' : addParams
  """
    output=\$(python3 ${projectDir}/bin/cluster_statistics.py ${flag} ${addParams} --toolName "${name}" "${newick}" "${sequences}" "${clusterFile}")
    echo ${name},\$output > ${name}_stats.out
  """

}

process merge_evaluation {
  label 'evalMerger'
  publishDir "${params.output}/${params.eval_output}", mode: 'copy', pattern: "${sequences.baseName}_summary.csv"

  input:
  path(evaluations)
  path(sequences)

  output:
  path "${sequences.baseName}_summary.csv"

  script:
  """
  echo "Algorithm,Number_of_Sequences,Number_of_Cluster,smallest_cluster,largest_cluster,average_cluster_size,median_cluster_size,Average_distance_to_nearest_centroid,number_of_unclustered_sequences" > "${sequences.baseName}_summary.csv"
  cat ${evaluations} >> "${sequences.baseName}_summary.csv"
  """
}