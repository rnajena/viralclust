/************************************************************************
* VSEARCH
*
* Cluster sequences with the clustering module of VSEARCH.
************************************************************************/

process vclust {
  label 'vclust'
  publishDir "${params.output}/${params.vsearch_output}", mode: 'copy', pattern: '*vsearch*'

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_vsearch.fasta"
    path "${sequences.baseName}_vsearch_cluster.uc"
    path "${sequences.baseName}_vsearch_cluster.uc.clstr"

  script:
  """
  vsearch ${addParams} --threads ${task.cpus} --cluster_fast ${sequences} --centroids ${sequences.baseName}_vsearch.fasta --uc ${sequences.baseName}_vsearch_cluster.uc
  python3 ${projectDir}/bin/vclust2cdhit.py ${sequences.baseName}_vsearch_cluster.uc
  """


}