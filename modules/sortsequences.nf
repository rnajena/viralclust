/************************************************************************
* SORTSEQUENCES
*
* Sort all input sequences (fasta) into positive
* strands (highest ORF density)
************************************************************************/

process sort_sequences {
  label 'sortseq'
  publishDir  "${params.output}/${params.sort_output}", mode: 'copy', pattern: '*_positive.fasta'

  input:
  path(sequences)
  val sort_off

  output:
  path "${sequences.baseName}_positive.fasta", emit: sort_result

  script:
  """
  python3 ${projectDir}/bin/sort_sequences.py "${sequences}" "${sort_off}" > "${sequences.baseName}_positive.fasta"
  """
}