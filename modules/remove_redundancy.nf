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


  script:
  """
    mmseqs easy-linclust --min-seq-id 1.0 "${sequences}" "${sequences.baseName}_nr" tmp
    mv "${sequences.baseName}_nr_rep_seq.fasta" "${sequences.baseName}_nr.fasta"
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