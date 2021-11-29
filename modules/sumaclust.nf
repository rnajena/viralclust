/************************************************************************
* SUMACLUST
*
* Cluster sequences with sumaclust.
************************************************************************/

process sumaclust {
  label 'sumaclust'
  publishDir "${params.output}/${params.sumaclust_output}", mode: 'copy', pattern: '*.fasta*'
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_sumaclust.fasta'


  input:
    path(sequences)
    val(addParams)
    val(goi)

  output:
    tuple val ("${params.output}/${params.sumaclust_output}"), path ("${sequences.baseName}_sumaclust.fasta"), path("${sequences.baseName}_sumaclust.fasta.clstr")
    // path "${sequences.baseName}_sumaclust.fasta.clstr", emit: sumaclust_cluster
    path "${sequences.baseName}_sumaclust_UNCLUSTERED.fasta"


  script:
  def GOI = goi != 'NO FILE' ? "${goi}" : ''
  """
    sumaclust "${addParams}" -p "${task.cpus}" "${sequences}"  > "${sequences.baseName}_sumaclust.fasta"

    python3 ${projectDir}/bin/suma2cdhit.py "${sequences.baseName}_sumaclust.fasta" ${GOI}

    awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < "${sequences.baseName}_sumaclust.fasta" | tail -n +2 > tmp.fasta
    grep 'cluster_center=True' tmp.fasta | grep -v 'cluster_weight=1;' | xargs -n1 -I% grep -A1 "%" tmp.fasta  > "${sequences.baseName}_sumaclust.fasta"
    grep 'cluster_center=True' tmp.fasta | grep 'cluster_weight=1;' | xargs -n1 -I% grep -A1 "%" tmp.fasta > "${sequences.baseName}_sumaclust_UNCLUSTERED.fasta"
    rm tmp.fasta

    if [ "{$GOI}" != 'NO FILE' ]; then
      for ID in \$(grep '>' ${GOI}); do
        grep -m 1 "\$ID" "${sequences.baseName}_sumaclust.fasta" || grep -A1 "\$ID" ${GOI}   >> "${sequences.baseName}_sumaclust.fasta"
      done
    fi

  """


}