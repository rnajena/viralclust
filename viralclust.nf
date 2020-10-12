#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
* Clustering of viral genomes based on different algorithms and metrices
*
* Author: kevin.lamkiewicz@uni-jena.de
*/

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $workflow.start"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
println "Configuration files:"
println "  $workflow.configFiles\u001B[0m"
println " "

if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

if ( params.profile ) {
  exit 1, "ERROR: --profile is WRONG use -profile"
}

if ( params.fasta == '' ) {
  exit 1, "ERROR: --fasta is a required parameter"
}

sortSequence = Channel.fromPath( workflow.projectDir + '/bin/sort_sequences.py', checkIfExists: true )
umap_hdbscan_script = Channel.fromPath( workflow.projectDir + '/bin/hdbscan_virus.py', checkIfExists: true )
umap_hdbscan_class = Channel.fromPath( workflow.projectDir + '/bin/ClusterViruses.py', checkIfExists: true )
sumaclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/suma2cdhit.py', checkIfExists: true )
cdhit2cdhit = Channel.fromPath( workflow.projectDir + '/bin/cdhit2goodcdhit.py', checkIfExists: true )
vclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/vclust2cdhit.py', checkIfExists: true )
prepareCSS = Channel.fromPath( workflow.projectDir + '/bin/prepare_css.py', checkIfExists: true )
clusterStats = Channel.fromPath( workflow.projectDir + '/bin/cluster_statistics.py', checkIfExists: true )

implicitTree = false
eval_params = ''
if (params.ncbi) {
  implicitTree = true
  eval_params = '--ncbi'
}

log.info """\
    VIRALCLUST -- CLUSTER YOUR VIRUSES
    ==================================
    Input File:             $params.fasta
    Output path:            $params.output
    CPUs used:              $params.cores
    ${msg -> if (params.tree | implicitTree) msg << "Tree will be calculated"}
    ${sw -> if (params.cdhit_params != '') sw << "cd-hit-est parameters:  ${params.cdhit_params}"}
    ${sw -> if (params.hdbscan_params != '') sw << "HDBscan parameters:     ${params.hdbscan_params}"}
    ${sw -> if (params.sumaclust_params != '') sw << "sumaclust parameters:     ${params.sumaclust_params}"}
    ${sw -> if (params.vsearch_params != '') sw << "vsearch parameters:     ${params.vsearch_params}"}


    """
    .stripIndent()

sequences = Channel.fromPath(params.fasta)

include { sort_sequences } from './modules/sortsequences'
include { remove_redundancy; cdhit } from './modules/cdhit'
include { hdbscan } from './modules/hdbscan'
include { sumaclust } from './modules/sumaclust'
include { vclust } from './modules/vsearch'
include { reverseComp } from './modules/reverseComp'

if (params.tree | implicitTree) {
  include { mafft } from './modules/mafft'
  include { fasttree } from './modules/fasttree'
  include { nwdisplay } from './modules/nwutils'
  include { prepare_css } from './modules/prepare_css'
  include { evaluate_cluster; merge_evaluation } from './modules/evaluate'
}


workflow {
  sort_sequences(sequences)
  remove_redundancy(sort_sequences.out.sort_result)

  hdbscan(remove_redundancy.out.nr_result, params.hdbscan_params)
  cdhit(remove_redundancy.out.nr_result, params.cdhit_params)
  sumaclust(remove_redundancy.out.nr_result, params.sumaclust_params)
  vclust(remove_redundancy.out.nr_result, params.vsearch_params)

  revCompChannel = hdbscan.out.hdbscan_result.concat(cdhit.out.cdhit_result, sumaclust.out.sumaclust_result, vclust.out.vclust_result)
  reverseComp(revCompChannel)

  if (params.tree | implicitTree) {
    mafft(revCompChannel)
    fasttree(mafft.out.mafft_result)
    nwdisplay(fasttree.out.fasttree_result)
    // mafft(remove_redundancy.out.nr_result)
    // fasttree(mafft.out.mafft_result)
    // colorChannel = vclust.out.vclust_cluster.concat(sumaclust.out.sumaclust_cluster, cdhit.out.cdhit_cluster, hdbscan.out.hdbscan_cluster)
    // prepare_css(colorChannel)
    // nwChannel = fasttree.out.fasttree_result.combine(prepare_css.out.css_cluster)
    // nwdisplay(nwChannel)

    // hdbEval = Channel.value('HDBscan').combine(hdbscan.out.hdbscan_cluster)
    // cdhitEval = Channel.value('cd-hit-est').combine(cdhit.out.cdhit_cluster)
    // sumaEval = Channel.value('sumaclust').combine(sumaclust.out.sumaclust_cluster)
    // vclustEval = Channel.value('vclust').combine(vclust.out.vclust_cluster)

    // clusterEval = hdbEval.concat(cdhitEval, sumaEval, vclustEval)
    // evalChannel = clusterEval.combine(fasttree.out.fasttree_result).combine(remove_redundancy.out.nr_result).combine(Channel.value(eval_params))


    // evaluate_cluster(evalChannel)
    // merge_evaluation(evaluate_cluster.out.eval_result.collect(), sequences)

  }


}