/************************************************************************
* CDHIT
*
* Cluster sequences by similarity -- here, remove 100% redundancy
************************************************************************/

process cdhit {
  label 'cdhit'
  publishDir "${params.output}/${params.cdhit_output}", mode: 'copy', pattern: "*_nr.fasta"

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_nr.fasta", emit: cdhit_result
    path "${sequences.baseName}_nr.error.log"

  script:  
  """
    cd-hit-est ${addParams} -i ${sequences} -o "${sequences.baseName}_nr.fasta"  2>"${sequences.baseName}_nr.error.log"
  """
}