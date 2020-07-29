/************************************************************************
* SUMACLUST
*
* Cluster sequences with sumaclust.
************************************************************************/

process sumaclust {
  label 'sumaclust'
  publishDir "${params.output}/${params.sumaclust_output}", mode: 'copy', pattern: '*sumaclust.fasta'
  

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_sumaclust.fasta", emit: sumaclust_result

  script:
  """
  sumaclust ${addParams} ${sequences} | grep "cluster_center=True" > "${sequences.baseName}_sumaclust.fasta"
  """


}