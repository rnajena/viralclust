/************************************************************************
* REMOVE REDUNDANCY
*
* Cluster 100% identical sequences to remove redundancy
************************************************************************/

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
    cd-hit-est -M 4000 -c 1 -i ${sequences} -o "${sequences.baseName}_nr.fasta"  2>"${sequences.baseName}_nr.error.log"
    sed -E "/>/ s/[:,()' ;]/_/g" "${sequences.baseName}_nr.fasta" > ${sequences.baseName}_renamed.fasta
    mv ${sequences.baseName}_renamed.fasta ${sequences.baseName}_nr.fasta
  """
}

process concat_goi {
  label 'concat_goi'
  publishDir "${params.output}/${params.nr_output}", mode: 'copy', overwrite: true, pattern: "*_nr.fasta"

  input:
    path(sequences)
    path(goi)

  output:
    path "${sequences.baseName}.fasta", emit: nr_result

  script:
  """
    cat "${goi}" "${sequences}" > tmp.fa 
    mv tmp.fa "${sequences}"
  """
}