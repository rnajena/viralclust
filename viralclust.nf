#!/usr/bin/env nextflow

nextflow.preview.dsl=2


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
println "Permanent cache directory:"
println "  $params.permanentCacheDir"
println "Configuration files:"
println "  $workflow.configFiles\u001B[0m"
println " "

if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPS to use: $params.max_cores\u001B[0m"
    println " "
}

if ( params.profile ) { 
  exit 1, "ERROR: --profile is WRONG use -profile" 
}

if ( params.fasta == '' ) {
  exit 1, "ERROR: --fasta is a required parameter"
}


umap_hdbscan_script = Channel.fromPath( workflow.projectDir + '/bin/viralClust.py', checkIfExists: true )


log.info """\
    VIRALCLUST -- CLUSTER YOUR VIRUSES
    ==================================
    Input File:             $params.fasta
    Output path:            $params.output
    CPUs used:              $params.cores

    cd-hit-est parameters:  $params.cdhit_params
    HDBscan parameters:     $params.hdbscan_params
    RAxML-ng parameters:    $params.raxmlng_params


    """
    .stripIndent()

sequences = Channel.fromPath(params.fasta)
// Channel.fromPath(params.fasta).set{sequences}

include { remove_redundancy; cdhit } from './modules/cdhit'
// include { rename_fasta_header } from './modules/rename_fasta_header'
//include { cdhit } from './modules/cdhit'
include { hdbscan } from './modules/hdbscan'
if (params.tree) {include { mafft } from './modules/mafft'}
if (params.tree) {include { raxmlng } from './modules/raxml-ng'}
//include { mafft } from './modules/mafft'
//include { raxmlng } from './modules/raxml-ng'

workflow {
  remove_redundancy(sequences)
  // rename_fasta_header(cdhit.out.cdhit_result)
  if (params.tree) {
    mafft(remove_redundancy.out.nr_result)
    raxmlng(mafft.out.mafft_result, params.raxmlng_params)
  }
  hdbscan(remove_redundancy.out.nr_result, params.hdbscan_params)
  cdhit(remove_redundancy.out.nr_result, params.cdhit_params)
  // uclust(cdhit.out.cdhit_result)
  // linclust(cdhit.out.cdhit_result)
  // applyMetric(hdb.out, uclust.out, linclust.out)

}