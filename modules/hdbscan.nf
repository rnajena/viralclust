/************************************************************************
* HDBSCAN
*
* Utilizing UMAP and HDBSCAN to cluster sequences based on their 
* k-mer vector representation and the cosine distance
************************************************************************/

process hdbscan {
  label 'hdbscan'
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*_hdbscan.fasta"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*.txt"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "cluster*.fasta"

  input:
    path(sequences)
    val(addParams)

  output:
    path "${sequences.baseName}_hdbscan.fasta", emit: hdbscan_result
    path "*"

  script:
  """
    python3 ${baseDir}/bin/viralClust.py -p ${task.cpus} ${addParams} ${sequences}
  """

}