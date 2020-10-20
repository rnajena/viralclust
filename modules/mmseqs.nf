/************************************************************************
* MMSEQS
*
* Use MMSEQs easy-linclust to cluster the input data
************************************************************************/

process mmseqs{
  label 'mmseqs'
  publishDir "${params.output}/${params.mmseqs_output}", mode: 'copy', pattern: "*_mmseqs*"

  input:
    path(sequences)
    val(addParams)

  output:
    tuple val("${params.output}/${params.mmseqs_output}"), path("${sequences.baseName}_mmseqs.fasta"), emit: mmseqs_result
    path "${sequences.baseName}_mmseqs.fasta.clstr", emit: mmseqs_cluster

  script:
  """
    mmseqs easy-linclust "${sequences}" "${sequences.baseName}_mmseqs" tmp
    mv ${sequences.baseName}_mmseqs_rep_seq.fasta ${sequences.baseName}_mmseqs.fasta

    python3 ${projectDir}/bin/mmseqs2cdhit.py ${sequences.baseName}_mmseqs_cluster.tsv "${sequences}"
    mv ${sequences.baseName}_mmseqs_cluster.tsv.clstr "${sequences.baseName}_mmseqs.fasta.clstr"

  """
}