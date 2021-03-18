/************************************************************************
* MAFFT
*
* Build a multiple sequence alignment using mafft.
************************************************************************/

process mafft {
  label 'mafft'
  publishDir "${params.output}/${params.mafft_output}", mode: 'copy', pattern: '*aln'

  input:
    tuple val(name), val(path), path(sequences), path(cluster)

  output:
    tuple val(name), path("${sequences.baseName}_mafft.aln"), emit: mafft_result

  script:
  """
  mafft --thread ${task.cpus} --reorder ${sequences} > "${sequences.baseName}_mafft.aln"
  """


}