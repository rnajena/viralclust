/************************************************************************
* RAXML-NG
*
* Inferring a ML tree based on a multiple sequence alignment
************************************************************************/

process raxmlng {
  label 'raxmlng'
  publishDir "${params.output}/${params.raxml-ng_output}", mode: 'copy'

  input:
    path(alignment)
    val(addParams)

  output:
    

  script:
  """
    raxml-ng --all --msa ${alignment} ${addParams} --threads ${params.cores}
  """


}