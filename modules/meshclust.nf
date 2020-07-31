/************************************************************************
* MESHCLUST
*
* Cluster sequences with Meshclust / Meshclust2
************************************************************************/

process meshclust {
  label 'meshclust'
  publishDir "${params.output}/${params.meshclust_output}", mode: 'copy', pattern: ''

  input:
    path(sequences)

  output:
    

  script:
  """
  meshclust ${sequences}
  """

}

process meshclust2 {
  label 'meshclust2'
  publishDir "${params.output}/${params.meshclust_output}", mode: 'copy', pattern: ''

  input:
    path(sequences)

  output:
    

  script:
  """
  meshclust2 ${sequences}
  """

}