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

if ( params.fasta == ''  && ! params.update_ncbi) {
  exit 1, "ERROR: --fasta or --update_ncbi is a required parameter.\nMake sure at least one is set."
}

sortSequence = Channel.fromPath( workflow.projectDir + '/bin/sort_sequences.py', checkIfExists: true )
rc_script = Channel.fromPath( workflow.projectDir + '/bin/reverse_complement.py', checkIfExists: true )
utils = Channel.fromPath( workflow.projectDir + '/bin/utils.py', checkIfExists: true )
ncbi_script = Channel.fromPath( workflow.projectDir + '/bin/get_ncbi_information.py', checkIfExists: true)
umap_hdbscan_script = Channel.fromPath( workflow.projectDir + '/bin/hdbscan_virus.py', checkIfExists: true )
umap_hdbscan_class = Channel.fromPath( workflow.projectDir + '/bin/ClusterViruses.py', checkIfExists: true )
sumaclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/suma2cdhit.py', checkIfExists: true )
cdhit2cdhit = Channel.fromPath( workflow.projectDir + '/bin/cdhit2goodcdhit.py', checkIfExists: true )
vclust2cdhit = Channel.fromPath( workflow.projectDir + '/bin/vclust2cdhit.py', checkIfExists: true )
mmseq2cdhit = Channel.fromPath( workflow.projectDir + '/bin/mmseqs2cdhit.py', checkIfExists: true )
clusterStats = Channel.fromPath( workflow.projectDir + '/bin/cluster_statistics.py', checkIfExists: true )
// ncbiMeta = Channel.fromPath ( workflow.projectDir + '/data/ncbi_metainfo.pkl', checkIfExists: true )

implicitEval = false
eval_params = ''
ncbiEval = Channel.from(false).combine(Channel.from(false))

if (params.ncbi) {
  implicitEval = true
  eval_params = "--ncbi"// ${projectDir}/data/ncbi_metainfo.pkl"
}

log.info """\
    VIRALCLUST -- CLUSTER YOUR VIRUSES
    ==================================
    Input File:             $params.fasta
    Output path:            $params.output
    CPUs used:              $params.cores
    ${msg -> if (params.eval | implicitEval) msg << "Tree will be calculated"}
    ${sw -> if (params.cdhit_params != '') sw << "cd-hit-est parameters:  ${params.cdhit_params}"}
    ${sw -> if (params.hdbscan_params != '') sw << "HDBscan parameters:     ${params.hdbscan_params}"}
    ${sw -> if (params.sumaclust_params != '') sw << "sumaclust parameters:     ${params.sumaclust_params}"}
    ${sw -> if (params.vclust_params != '') sw << "vclust parameters:     ${params.vclust_params}"}


    """
    .stripIndent()

if (params.fasta) {

  sequences = Channel.fromPath(params.fasta)

  include { sort_sequences } from './modules/sortsequences'
  include { remove_redundancy; cdhit } from './modules/cdhit'
  include { hdbscan } from './modules/hdbscan'
  include { sumaclust } from './modules/sumaclust'
  include { vclust } from './modules/vclust'
  include { mmseqs } from './modules/mmseqs'
  include { reverseComp } from './modules/reverseComp'

  if (params.eval | implicitEval) {
    include { mafft } from './modules/mafft'
    include { fasttree } from './modules/fasttree'
    include { nwdisplay } from './modules/nwutils'
    include { evaluate_cluster; merge_evaluation } from './modules/evaluate'
  }

  if (params.ncbi) {
    include { get_ncbi_meta } from './modules/ncbi_meta'
    ncbi_metainfo_ch = file("${workflow.projectDir}/data/ncbi_metainfo.pkl")
    if (! ncbi_metainfo_ch.exists() & ! params.update_ncbi) {
      include { update_ncbi_metainfo } from './modules/update_ncbi_metainfo'

      log.warn """ \
      No NCBI meta information database found.
      Database will be downloaded and stored for this and future runs!
      """.stripIndent()
    } else {
      ncbiMeta = Channel.fromPath ( workflow.projectDir + '/data/ncbi_metainfo.pkl')
    }
  }
}


if (params.update_ncbi) {
  include {update_ncbi_metainfo} from './modules/update_ncbi_metainfo'
}

workflow update_metadata {
  update_ncbi_metainfo()
}

workflow annotate_metadata {
  take:
    non_redundant_ch

  main:
    get_ncbi_meta(non_redundant_ch)
    ncbiEval = Channel.from(eval_params).combine(ncbiMeta)

  emit:
    ncbiEval
}

workflow preprocessing {
  main:
    sort_sequences(sequences)
    remove_redundancy(sort_sequences.out.sort_result)
    non_redundant_ch = remove_redundancy.out.nr_result
  emit:
    non_redundant_ch
}

workflow clustering {

  take:
    non_redundant_ch

  main:
    // if (params.ncbi) {
    //   annotate_metadata(non_redundant_ch)
    //   ncbiEval = annotate_metadata.out.ncbiEval
    // }

    hdbscan(non_redundant_ch, params.hdbscan_params)
    cdhit(non_redundant_ch, params.cdhit_params)
    sumaclust(non_redundant_ch, params.sumaclust_params)
    vclust(non_redundant_ch, params.vclust_params)
    mmseqs(non_redundant_ch, params.vclust_params)


    hdbRC = Channel.value('HDBSCAN').combine(hdbscan.out.hdbscan_result)
    cdhitRC = Channel.value('cd-hit-est').combine(cdhit.out.cdhit_result)
    sumaRC = Channel.value('sumaclust').combine(sumaclust.out.sumaclust_result)
    vclustRC = Channel.value('vclust').combine(vclust.out.vclust_result)
    mmseqsRC = Channel.value('MMseqs2').combine(mmseqs.out.mmseqs_result)
    revCompChannel = hdbRC.concat(cdhitRC, sumaRC, vclustRC, mmseqsRC)
    reverseComp(revCompChannel)

    hdbEval = Channel.value('HDBSCAN').combine(hdbscan.out.hdbscan_cluster)
    cdhitEval = Channel.value('cd-hit-est').combine(cdhit.out.cdhit_cluster)
    sumaEval = Channel.value('sumaclust').combine(sumaclust.out.sumaclust_cluster)
    vclustEval = Channel.value('vclust').combine(vclust.out.vclust_cluster)
    mmseqsEval = Channel.value('MMseqs2').combine(mmseqs.out.mmseqs_cluster)
    clusterEval = hdbEval.concat(cdhitEval, sumaEval, vclustEval, mmseqsEval)

  emit:
    revCompChannel
    clusterEval
}

workflow evaluation {
  take:
    revCompChannel
    clusterEval
    non_redundant_ch

  main:
    mafft(revCompChannel)
    fasttree(mafft.out.mafft_result)
    nwdisplay(fasttree.out.fasttree_result)

    if (params.ncbi) {
      if (! ncbi_metainfo_ch.exists() & ! params.update_ncbi) {
        update_metadata()
        ncbiMeta = Channel.fromPath ( workflow.projectDir + '/data/ncbi_metainfo.pkl')
      }
      ncbiEval = Channel.from(eval_params).combine(ncbiMeta)
    }

    evalChannel = clusterEval.join(fasttree.out.fasttree_result).combine(non_redundant_ch).combine(ncbiEval)
    evaluate_cluster(evalChannel)
    merge_evaluation(evaluate_cluster.out.eval_result.collect(), sequences)

    //sequences = Channel.fromPath(params.fasta)
    evaluate_cluster.out.warning.first().subscribe{

      log.warn """\

      ##########################################################
      NCBI meta information is older than 90 days.
      Please consider updating using the following command:
        nextflow run viralclust.nf --update_ncbi
      ##########################################################
      """.stripIndent()
    }
}

workflow {

  if (params.update_ncbi) {
    update_metadata()
  }

  if (params.fasta) {
    preprocessing()
    clustering(preprocessing.out.non_redundant_ch)
  }

  if (params.eval | implicitEval) {
    evaluation(clustering.out.revCompChannel, clustering.out.clusterEval, preprocessing.out.non_redundant_ch)
  }
  // sort_sequences(sequences)
  // remove_redundancy(sort_sequences.out.sort_result)

  // if (params.ncbi) {
  //   get_ncbi_meta(remove_redundancy.out.nr_result)
  //   ncbiEval = Channel.from(eval_params).combine(get_ncbi_meta.out.pkl_ncbi)
  // }

  // hdbscan(remove_redundancy.out.nr_result, params.hdbscan_params)
  // cdhit(remove_redundancy.out.nr_result, params.cdhit_params)
  // sumaclust(remove_redundancy.out.nr_result, params.sumaclust_params)
  // vclust(remove_redundancy.out.nr_result, params.vclust_params)
  // mmseqs(remove_redundancy.out.nr_result, params.vclust_params)


  // hdbRC = Channel.value('HDBSCAN').combine(hdbscan.out.hdbscan_result)
  // cdhitRC = Channel.value('cd-hit-est').combine(cdhit.out.cdhit_result)
  // sumaRC = Channel.value('sumaclust').combine(sumaclust.out.sumaclust_result)
  // vclustRC = Channel.value('vclust').combine(vclust.out.vclust_result)
  // mmseqsRC = Channel.value('MMseqs2').combine(mmseqs.out.mmseqs_result)
  // revCompChannel = hdbRC.concat(cdhitRC, sumaRC, vclustRC, mmseqsRC)
  // reverseComp(revCompChannel)

  // if (params.eval | implicitEval) {
  //   mafft(revCompChannel)
  //   fasttree(mafft.out.mafft_result)
  //   nwdisplay(fasttree.out.fasttree_result)

  //   hdbEval = Channel.value('HDBSCAN').combine(hdbscan.out.hdbscan_cluster)
  //   cdhitEval = Channel.value('cd-hit-est').combine(cdhit.out.cdhit_cluster)
  //   sumaEval = Channel.value('sumaclust').combine(sumaclust.out.sumaclust_cluster)
  //   vclustEval = Channel.value('vclust').combine(vclust.out.vclust_cluster)
  //   mmseqsEval = Channel.value('MMseqs2').combine(mmseqs.out.mmseqs_cluster)

  //   clusterEval = hdbEval.concat(cdhitEval, sumaEval, vclustEval, mmseqsEval)

  //   evalChannel = clusterEval.join(fasttree.out.fasttree_result).combine(remove_redundancy.out.nr_result).combine(ncbiEval)
  //   evaluate_cluster(evalChannel)
  //   merge_evaluation(evaluate_cluster.out.eval_result.collect(), sequences)
  //}


}