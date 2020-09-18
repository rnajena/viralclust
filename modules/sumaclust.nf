/************************************************************************
* SUMACLUST
*
* Cluster sequences with sumaclust.
************************************************************************/

process sumaclust {
  label 'sumaclust'
  publishDir "${params.output}/${params.sumaclust_output}", mode: 'copy', pattern: '*sumaclust.fasta*'
  

  input:
    path(sequences)
    val(addParams)

  output:
    tuple val ("${params.output}/${params.sumaclust_output}"), path ("${sequences.baseName}_sumaclust.fasta"), emit: sumaclust_result
    path "${sequences.baseName}_sumaclust.fasta.clstr", emit: sumaclust_cluster
    

  script:
  """
  sumaclust ${addParams} ${sequences}  > "${sequences.baseName}_sumaclust.fasta"
  python3 ${projectDir}/bin/suma2cdhit.py "${sequences.baseName}_sumaclust.fasta"
  
  """


}