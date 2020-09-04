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
umap_hdbscan_script = Channel.fromPath( workflow.projectDir + '/bin/viralClust.py', checkIfExists: true )
umap_hdbscan_class = Channel.fromPath( workflow.projectDir + '/bin/ClusterViruses.py', checkIfExists: true )
sumaclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/suma2cdhit.py', checkIfExists: true )
vclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/vclust2cdhit.py', checkIfExists: true )
prepareCSS = Channel.fromPath( workflow.projectDir + '/bin/prepare_css.py', checkIfExists: true )

log.info """\
    VIRALCLUST -- CLUSTER YOUR VIRUSES
    ==================================
    Input File:             $params.fasta
    Output path:            $params.output
    CPUs used:              $params.cores
    ${msg -> if (params.tree) msg << "Tree will be calculated"}
    ${sw -> if (params.cdhit_params != '') sw << "cd-hit-est parameters:  ${params.cdhit_params}"}
    ${sw -> if (params.hdbscan_params != '') sw << "HDBscan parameters:     ${params.hdbscan_params}"}
    ${sw -> if (params.sumaclust_params != '') sw << "sumaclust parameters:     ${params.sumaclust_params}"}
    ${sw -> if (params.vsearch_params != '') sw << "vsearch parameters:     ${params.vsearch_params}"}
    

    """
    .stripIndent()
// ${sw -> if (params.tree) sw << "RAxML-ng parameters:    ${params.raxmlng_params}"}

sequences = Channel.fromPath(params.fasta)
// Channel.fromPath(params.fasta).set{sequences}

include { sort_sequences } from './modules/sortsequences'
include { remove_redundancy; cdhit } from './modules/cdhit'
include { hdbscan } from './modules/hdbscan'
include { sumaclust } from './modules/sumaclust'
include { vclust } from './modules/vsearch'

if (params.tree) {
  include { mafft } from './modules/mafft'
  include { fasttree } from './modules/fasttree'
  include { raxmlng } from './modules/raxml-ng'
  include { nwdisplay } from './modules/nwutils'
  include { prepare_css } from './modules/prepare_css'
}


workflow {
  sort_sequences(sequences)
  remove_redundancy(sort_sequences.out.sort_result)

  hdbscan(remove_redundancy.out.nr_result, params.hdbscan_params)
  cdhit(remove_redundancy.out.nr_result, params.cdhit_params)
  sumaclust(remove_redundancy.out.nr_result, params.sumaclust_params)
  vclust(remove_redundancy.out.nr_result, params.vsearch_params)


  if (params.tree) {
    mafft(remove_redundancy.out.nr_result)
    fasttree(mafft.out.mafft_result)
    colorChannel = vclust.out.vclust_cluster.concat(sumaclust.out.sumaclust_cluster, cdhit.out.cdhit_cluster, hdbscan.out.hdbscan_cluster)
    //colorChannel.view()
    prepare_css(colorChannel)
    //prepare_css.out.css_cluster.collate(3).view()
    nwChannel = fasttree.out.fasttree_result.combine(prepare_css.out.css_cluster)
    nwdisplay(nwChannel)
    //nwChannel.view()
    //nwdisplay(fasttree.out.fasttree_result, prepare_css.out.css_cluster.collect())
    //prepare_css.out.css_cluster.subscribe( nwdisplay(fasttree.out.fasttree_result, prepare_css.out.css_cluster, prepare_css.out.css_ornaments) )


    // prepare_css.out.css_cluster.view()
    // prepare_css.out.css_ornaments.view()
    // nwdisplay(fasttree.out.fasttree_result, prepare_css.out.css_cluster)
    //raxmlng(mafft.out.mafft_result, params.raxmlng_params)
  }


}