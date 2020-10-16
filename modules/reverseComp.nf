/************************************************************************
* Reverse Comp
*
* Generate reverse complementary sequences based on input records.
************************************************************************/

process reverseComp {
  label 'revcomp'
  publishDir "${path}", mode: 'copy', pattern: '*.fasta'


  input:
    tuple val(name), val(path), path(sequences)


  output:
    tuple val(name), path(outpath)

  script:
  outpath = "${sequences}".reverse().replaceFirst("positive".reverse(),"negative".reverse()).reverse()
  """
  python3 ${projectDir}/bin/reverse_complement.py "${sequences}"

  """


}