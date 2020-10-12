/************************************************************************
* NEWICK_UTILITIES
*
* Manipulate and visualize NEWICK trees
************************************************************************/

process nwdisplay {
  label 'nwdisplay'
  publishDir "${params.output}/${params.nwdisplay_output}", mode: 'copy', pattern: "*_nwdisplay*"
  publishDir "${params.output}/${params.eval_output}", mode: 'copy', pattern: "*pdf"

  input:
    //tuple path(newick), val(name), path(css), path(ornaments)
    path(newick)

  output:
    //path "${name}_nwdisplay.svg", emit: nwdisplay_result
    //path "${name}_nwdisplay.pdf"
    path "${newick}_nwdisplay.svg", emit: nwdisplay_result
    path "${newick}_nwdisplay.pdf"
  

  script:
  """
  nw_display -v 25 -i "font-size:2" -l "font-size:6;font-family:helvetica;font-style:italic" -Il -w 5000 -r -b "opacity:0" -s ${newick} > "${newick}_nwdisplay.svg"
  rsvg-convert -f pdf -o "${newick}_nwdisplay.pdf" ${newick}_nwdisplay.svg

  """
}