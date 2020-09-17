/************************************************************************
* SUMACLUST
*
* Cluster sequences with sumaclust.
************************************************************************/

process reverseComp {
  label 'revcomp'
  publishDir "${path}", mode: 'copy', pattern: '*.fasta'
  

  input:
    tuple val(path), path(sequences)
    

  output:
    path(outpath)  

  script:
  outpath = "${sequences}".reverse().replaceFirst("positive".reverse(),"negative".reverse()).reverse()
  """
  python3 ${projectDir}/bin/reverse_complement.py "${sequences}"

  """


}