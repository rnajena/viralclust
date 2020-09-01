/************************************************************************
* CD-HIT-EST
*
* Cluster sequences by similarity
************************************************************************/

process cdhit {
  label 'cdhit'
  publishDir "${params.output}/${params.cdhit_output}", mode: 'copy', pattern: '*_cdhitest.fasta*'

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_cdhitest.fasta", emit: cdhit_result
    path "${sequences.baseName}_cdhitest.fasta.clstr"

  script:  
  """
    cd-hit-est ${addParams} -i ${sequences} -o "${sequences.baseName}_cdhitest.fasta"
  """
}


process remove_redundancy {
  label 'remove'
  publishDir "${params.output}/${params.nr_output}", mode: 'copy', pattern: "*_nr.fasta"

  input:
    path(sequences)

  output:
    path "${sequences.baseName}_nr.fasta", emit: nr_result
    path "${sequences.baseName}_nr.error.log"

  script:  
  """
    cd-hit-est -c 1 -i ${sequences} -o "${sequences.baseName}_nr.fasta"  2>"${sequences.baseName}_nr.error.log"
    sed -E "/>/ s/[:,()' ;]/_/g" "${sequences.baseName}_nr.fasta" > ${sequences.baseName}_renamed.fasta
    mv ${sequences.baseName}_renamed.fasta ${sequences.baseName}_nr.fasta
  """
}