/************************************************************************
* CD-HIT-EST
*
* Cluster sequences by similarity
************************************************************************/

process cdhit {
  label 'cdhit'
  publishDir "${params.output}/${params.cdhit_output}", mode: 'copy', pattern: '*_cdhitest.fasta*'
  publishDir "${params.output}/${params.cdhit_output}", mode: 'copy', pattern: '*UNCLUSTERED*'
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_cdhitest.fasta'

  input:
    path(sequences)
    val(addParams)
    val(goi)

  output:
    tuple val("${params.output}/${params.cdhit_output}"), path("${sequences.baseName}_cdhitest.fasta"), path("${sequences.baseName}_cdhitest.fasta.clstr")
    //path "${sequences.baseName}_cdhitest.fasta.clstr", emit: cdhit_cluster
    path "${sequences.baseName}_cdhitest_UNCLUSTERED.fasta"

  script:
  def GOI = goi != 'NO FILE' ? "${goi}" : ''
  """
    cd-hit-est ${addParams} -T "${task.cpus}" -i ${sequences} -o "${sequences.baseName}_cdhitest.fasta"

    python3 ${baseDir}/bin/cdhit2goodcdhit.py "${sequences.baseName}_cdhitest.fasta.clstr" ${sequences} ${GOI} > tmp.clstr
    mv tmp.clstr "${sequences.baseName}_cdhitest.fasta.clstr"
    python3 ${baseDir}/bin/filter_unclustered.py "${sequences.baseName}_cdhitest.fasta" "${sequences.baseName}_cdhitest.fasta.clstr"
    mv "${sequences.baseName}_cdhitest.fastaTEST" "${sequences.baseName}_cdhitest.fasta"

    if [ "{$GOI}" != 'NO FILE' ]; then
      for ID in \$(grep '>' ${GOI}); do
        grep -m 1 "\$ID" "${sequences.baseName}_cdhitest.fasta" || grep -A1 "\$ID" ${GOI}  >> "${sequences.baseName}_cdhitest.fasta"
      done 
    fi

  """
}
