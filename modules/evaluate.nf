/************************************************************************
* EVALUATE
*
* Calculate basic statistics for each cluster algorithm
************************************************************************/

process evaluate_cluster {
  label 'evaluate'

  input:
    tuple val(name), path(clusterFile), path(newick), path(sequences)

  output:
    path "${name}_stats.out", emit: eval_result

  script:
  """
    echo ${name}, \$(python3 ${projectDir}/bin/cluster_statistics.py "${newick}" "${sequences}" "${clusterFile}") >> ${name}_stats.out    
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
  echo "Algorithm, Number of Sequences, Number of Cluster, smallest cluster, largest cluster, average cluster size, median cluster size , Average distance to nearest centroid, number of unclustered sequences" > "${sequences.baseName}_summary.csv"
  cat ${evaluations} >> "${sequences.baseName}_summary.csv"
  """
}