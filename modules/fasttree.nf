/************************************************************************
* FASTTREE
*
* Build a phylogenetic tree with fasttree.
************************************************************************/

process fasttree {
  label 'fasttree'
  publishDir "${params.output}/${params.fasttree_output}", mode: 'copy', pattern: '*_fasttree.nwk'

  input:
    tuple val(name), path(alignment)

  output:
    tuple val(name), path("${alignment.baseName}_fasttree.nwk"), emit: fasttree_result

  script:
  """
  FastTreeMP -gtr -nt ${alignment} > "${alignment.baseName}_fasttree.nwk"
  """


}