/************************************************************************
* HDBSCAN
*
* Utilizing UMAP and HDBSCAN to cluster sequences based on their
* k-mer vector representation and the cosine distance
************************************************************************/

process hdbscan {
  label 'hdbscan'
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*_hdbscan.fasta"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*.clstr"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*log"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*_hdbscan_UNCLUSTERED.fasta"
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_hdbscan.fasta'


  input:
    path(sequences)
    val(addParams)
    val(goi)

  output:
    tuple val("${params.output}/${params.hdbscan_output}"), path("${sequences.baseName}_hdbscan.fasta"), emit: hdbscan_result
    path "${sequences.baseName}_hdbscan.fasta.clstr", emit: hdbscan_cluster
    path "${sequences.baseName}_hdbscan_UNCLUSTERED.fasta"
    path "hdbscan.log"

  script:
  def GOI = goi != 'NO FILE' ? "${goi}" : ''
  """
    python3 ${baseDir}/bin/hdbscan_virus.py -v -p ${task.cpus} ${addParams} ${sequences} ${GOI} 2> hdbscan.log
    mv cluster-1.fasta  "${sequences.baseName}_hdbscan_UNCLUSTERED.fasta"

  """

}