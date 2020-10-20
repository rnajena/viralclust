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
    tuple val(name), path(newick)

  output:
    path "${newick}_nwdisplay.svg", emit: nwdisplay_result
    path "${newick}_nwdisplay.pdf"

  script:
  """
  nw_reroot "${newick}" | nw_display -v 20 -i "font-size:3" -l "font-size:4;font-family:helvetica;font-style:italic" -b "font-size:3;opacity:0" -s - > "${newick}_nwdisplay.svg"
  rsvg-convert -f pdf -o "${newick}_nwdisplay.pdf" ${newick}_nwdisplay.svg

  """
}