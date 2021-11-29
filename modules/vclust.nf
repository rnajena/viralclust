/************************************************************************
* VCLUST
*
* Cluster sequences with the clustering module of VSEARCH.
************************************************************************/

process vclust {
  label 'vclust'
  publishDir "${params.output}/${params.vclust_output}", mode: 'copy', pattern: '*vclust*'
  publishDir "${params.output}/${params.vclust_output}", mode: 'copy', pattern: '*UNCLUSTERED*'
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_vclust.fasta'

  input:
    path(sequences)
    val(addParams)
    val(goi)

  output:
    tuple val("${params.output}/${params.vclust_output}"), path("${sequences.baseName}_vclust.fasta"), path("${sequences.baseName}_vclust_cluster.uc.clstr")
    path "${sequences.baseName}_vclust_cluster.uc"
    // path "${sequences.baseName}_vclust_cluster.uc.clstr", emit: vclust_cluster
    path "${sequences.baseName}_vclust_UNCLUSTERED.fasta"

  script:
  def GOI = goi != 'NO FILE' ? "${goi}" : ''
  """
    vsearch ${addParams} --threads ${task.cpus} --cluster_fast ${sequences} --centroids ${sequences.baseName}_vclust.fasta --uc ${sequences.baseName}_vclust_cluster.uc
    if [ "{$GOI}" != 'NO FILE' ]; then
      for ID in \$(grep '>' ${GOI}); do
        grep -m 1 "\$ID" "${sequences.baseName}_vclust.fasta" || grep -A1 "\$ID" ${GOI} >> "${sequences.baseName}_vclust.fasta"
      done
    fi


    python3 ${projectDir}/bin/vclust2cdhit.py ${sequences.baseName}_vclust_cluster.uc ${GOI}
    python3 ${projectDir}/bin/filter_unclustered.py "${sequences.baseName}_vclust.fasta" "${sequences.baseName}_vclust_cluster.uc.clstr"
    mv "${sequences.baseName}_vclust.fastaTEST" "${sequences.baseName}_vclust.fasta"

  """


}