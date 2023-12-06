/************************************************************************
* HTML
*
* Generate html page with the results.
************************************************************************/

process generate_html {
  label 'generate_html'
  publishDir "${params.output}/${params.summary_output}", mode: 'copy'

  input:
    path(clstrFile)

  output:
    path("web")

  script:
  """
    cp -r "${baseDir}/web/" .

    generate_html.py . web/

  """


}