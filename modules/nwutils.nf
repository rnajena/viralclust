/************************************************************************
* NEWICK_UTILITIES
*
* Manipulate and visualize NEWICK trees
************************************************************************/

process nwdisplay {
  label 'nwdisplay'
  publishDir "${params.output}/${params.nwdisplay_output}", mode: 'copy', pattern: "*_nwdisplay*"

  input:
    path(newick)

  output:
    path "${newick.baseName}_nwdisplay.svg", emit: nwdisplay_result
    path "${newick.baseName}_nwdisplay.pdf"

  script:
  """
  nw_display -v 25 -i "font-size:2" -l "font-size:6;font-family:helvetica;font-style:italic" -Il -w 5000 -r -b "opacity:0" -s ${newick} > "${newick.baseName}_nwdisplay.svg"
  rsvg-convert -f pdf -o "${newick.baseName}_nwdisplay.pdf" ${newick.baseName}_nwdisplay.svg

  """


}