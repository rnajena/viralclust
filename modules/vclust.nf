/************************************************************************
* VCLUST
*
* Cluster sequences with the clustering module of VSEARCH.
************************************************************************/

process vclust {
  label 'vclust'
  publishDir "${params.output}/${params.vclust_output}", mode: 'copy', pattern: '*vclust*'
  publishDir "${params.output}/${params.vclust_output}", mode: 'copy', pattern: '*UNCLUSTERED*'
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_vclust.fasta'

  input:
    path(sequences)
    val(addParams)

  output:
    tuple val("${params.output}/${params.vclust_output}"), path ("${sequences.baseName}_vclust.fasta"), emit: vclust_result
    path "${sequences.baseName}_vclust_cluster.uc"
    path "${sequences.baseName}_vclust_cluster.uc.clstr", emit: vclust_cluster
    path "${sequences.baseName}_vclust_UNCLUSTERED.fasta"

  script:
  """
  vsearch ${addParams} --threads ${task.cpus} --cluster_fast ${sequences} --centroids ${sequences.baseName}_vclust.fasta --uc ${sequences.baseName}_vclust_cluster.uc
  python3 ${projectDir}/bin/vclust2cdhit.py ${sequences.baseName}_vclust_cluster.uc
  python3 ${projectDir}/bin/filter_unclustered.py "${sequences.baseName}_vclust.fasta" "${sequences.baseName}_vclust_cluster.uc.clstr"

  """


}