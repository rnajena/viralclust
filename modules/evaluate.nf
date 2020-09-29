/************************************************************************
* EVALUATE
*
* Calculate basic statistics for each cluster algorithm
************************************************************************/

process evaluate_cluster {
  label 'evaluate'

  input:
    tuple val(name), path(newick), path(sequences), path(centroids), path(clusterFile)

  output:

  script:
  """
    python3 ${projectDir}/bin/cluster_statistics.py "${newick}" "${sequences}" "${centroids}" "${clusterFile}"
  """

}