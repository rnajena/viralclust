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
    tuple val("${params.output}/${params.vsearch_output}"), path ("${sequences.baseName}_vsearch.fasta"), emit: vclust_result
    path "${sequences.baseName}_vsearch_cluster.uc"
    path "${sequences.baseName}_vsearch_cluster.uc.clstr", emit: vclust_cluster

  script:
  """
  vsearch ${addParams} --threads ${task.cpus} --cluster_fast ${sequences} --centroids ${sequences.baseName}_vsearch.fasta --uc ${sequences.baseName}_vsearch_cluster.uc
  python3 ${projectDir}/bin/vclust2cdhit.py ${sequences.baseName}_vsearch_cluster.uc

  """


}