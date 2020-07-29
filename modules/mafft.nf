/************************************************************************
* MAFFT
*
* Build a multiple sequence alignment using mafft.
************************************************************************/

process mafft {
  label 'mafft'
  publishDir "${params.output}/${params.mafft_output}", mode: 'copy', pattern: '*aln'

  input:
    path(sequences)

  output:
    path "${sequences.baseName}_mafft.aln", emit: mafft_result

  script:
  """
  mafft --thread ${params.cores} --reorder ${sequences} > "${sequences.baseName}_mafft.aln"
  """


}