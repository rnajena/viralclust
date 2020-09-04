/************************************************************************
* NEWICK_UTILITIES
*
* Manipulate and visualize NEWICK trees
************************************************************************/

process nwdisplay {
  label 'nwdisplay'
  publishDir "${params.output}/${params.nwdisplay_output}", mode: 'copy', pattern: "*_nwdisplay*"

  input:
    tuple path(newick), val(name), path(css), path(ornaments)

  output:
    path "${name}_nwdisplay.svg", emit: nwdisplay_result
    path "${name}_nwdisplay.pdf"

  

  script:
  """
  nw_display -v 25 -i "font-size:2" -l "font-size:6;font-family:helvetica;font-style:italic" -Il -w 5000 -r -b "opacity:0" -s -c ${css} -o ${ornaments}  ${newick} > "${name}_nwdisplay.svg"
  rsvg-convert -f pdf -o "${name}_nwdisplay.pdf" ${name}_nwdisplay.svg

  """
}