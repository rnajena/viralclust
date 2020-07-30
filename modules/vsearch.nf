/************************************************************************
* VSEARCH
*
* Cluster sequences with the clustering module of VSEARCH.
************************************************************************/

process vsearch {
  label 'vsearch'
  publishDir "${params.output}/${params.vsearch_output}", mode: 'copy', pattern: '*vsearch*'

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_vsearch.fasta"
    path "${sequences.baseName}_vsearch_cluster.uc"

  script:
  """
  vsearch ${addParams} --cluster_fast ${sequences} --centroids ${sequences.baseName}_vsearch.fasta --uc ${sequences.baseName}_vsearch_cluster.uc
  """


}