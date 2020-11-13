#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
* Clustering of viral genomes based on different algorithms and metrices
*
* Author: kevin.lamkiewicz@uni-jena.de
*/

if (params.help) {
  exit 0, helpMSG()
}

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
//umap_hdbscan_class = Channel.fromPath( workflow.projectDir + '/bin/ClusterViruses.py', checkIfExists: true )
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

  if (params.goi) {
    goi = Channel.fromPath(params.goi)
    include { sort_sequences as goi_sorter } from './modules/sortsequences'
    include { concat_goi } from './modules/remove_redundancy'
  }

  include { sort_sequences } from './modules/sortsequences'
  include { remove_redundancy } from './modules/remove_redundancy'
  include { cdhit } from './modules/cdhit'
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
    ncbi_metainfo_ch = file("${workflow.projectDir}/${params.permanentCacheDir}/ncbi_metainfo.pkl")
    if (! ncbi_metainfo_ch.exists() & ! params.update_ncbi) {
      include { update_ncbi_metainfo } from './modules/update_ncbi_metainfo'

      log.warn """ \
      No NCBI meta information database found.
      Database will be downloaded and stored for this and future runs!
      """.stripIndent()
    } else {
      ncbiMeta = Channel.fromPath ( workflow.projectDir + params.permanentCacheDir + '/ncbi_metainfo.pkl')
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
    if (params.goi) {
      goi_sorter(goi)
      goiSorted = goi_sorter.out.sort_result
    } else {
      goiSorted = 'NO FILE'
    }

    remove_redundancy(sort_sequences.out.sort_result)
    non_redundant_ch = remove_redundancy.out.nr_result
    if (params.goi) {
      concat_goi(remove_redundancy.out.nr_result, goiSorted)
      non_redundant_ch = concat_goi.out.nr_result
    }
  emit:
    non_redundant_ch
    goiSorted
}

workflow clustering {

  take:
    non_redundant_ch
    goiSorted

  main:

    hdbscan(non_redundant_ch, params.hdbscan_params, goiSorted)
    cdhit(non_redundant_ch, params.cdhit_params, goiSorted)
    sumaclust(non_redundant_ch, params.sumaclust_params, goiSorted)
    vclust(non_redundant_ch, params.vclust_params, goiSorted)
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
        ncbiMeta = ncbiMeta = Channel.fromPath ( workflow.projectDir + params.permanentCacheDir + '/ncbi_metainfo.pkl')
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
    update_metadata(params.permanentCacheDir)
  }

  if (params.fasta) {
    preprocessing()
    clustering(preprocessing.out.non_redundant_ch, preprocessing.out.goiSorted)
  }

  if (params.eval | implicitEval) {
    evaluation(clustering.out.revCompChannel, clustering.out.clusterEval, preprocessing.out.non_redundant_ch)
  }
}

def helpMSG() {
    c_reset = "\033[0m";
    c_red = "\033[1;31m"
    c_green = "\033[1;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_green}Welcome to ViralClust - your pipeline to cluster viral genome sequences once and for all!${c_reset}
    ____________________________________________________________________________________________

    ${c_yellow}Usage example:${c_reset}
    nextflow run viralclust.nf --update_ncbi

    or

    nextflow run viralclust.nf --fasta "genomes.fasta"

    or both

    nextflow run viralclust.nf --update_ncbi --fasta "genomes.fasta"

    ____________________________________________________________________________________________

    ${c_yellow}Mandatory Input:${c_reset}
    ${c_green}--fasta PATH${c_reset}                      Path to a multiple fasta sequence file, storing all genomes that shall be clustered.
                                      Usually, this parameter has to be set, unless the parameter ${c_green}--ncbi_update${c_reset} has been set.

    ${c_yellow}Optional Input:${c_reset}
    ${c_green}--goi PATH${c_reset}                        Path to a (multiple) fasta sequence file with genomes that have to end
                                      up in the final set of representative genomes, e.g. strains of your lab that are
                                      of special interest. This parameter is optional.
    ____________________________________________________________________________________________

    ${c_yellow}Options:${c_reset}
    ${c_green}--eval${c_reset}                            After clustering, calculate basic statistics of clustering results. For each
                                      tool, the minimum, maximum, average and median cluster sizes are calculated,
                                      as well as the average distance of two representative genomes.

    ${c_green}--ncbi${c_reset}                            Additionally to the evaluation performed by ${c_green}--eval${c_reset}, NCBI metainformation
                                      is included for all genomes of the input set. Therefore, the identifier of fasta records are
                                      scanned for GenBank accession IDs, which are then used to retrieve information about the taxonomy,
                                      accession date and accession country of a sequence. Implicitly calls ${c_green}--eval${c_reset}.
                                      ${c_red}Attention:${c_reset} If no database is available at ${params.permanentCacheDir}, setting this flag
                                      implicitly sets ${c_green}--ncbi_update${c_reset}.

    ${c_green}--ncbi_update${c_reset}                     Downloads all current GenBank entries from the NCBI FTP server and processes the data to
                                      the databank stored at ${params.permanentCacheDir}.

    ${c_yellow}Cluster options:${c_reset}
    ${c_green}--cdhit_params${c_reset}                    Additional parameters for CD-HIT-EST cluster analysis. [default $params.cdhit_params]
                                      For more information and options, we refer to the CD-HIT manual.

    ${c_green}--hdbscan_params${c_reset}                  Additional parameters for HDBscan cluster analysis. [default $params.hdbscan_params]
                                      For more information and options, please use
                                      ${c_green}nextflow run viralclust.nf --hdbscan_help${c_reset} or ${c_green}python3 bin/hdbscan_virus.py -h${c_reset}.

    ${c_green}--sumaclust_params${c_reset}                Additional parameters for sumaclust cluster analysis. [default $params.sumaclust_params]
                                      For more information and options, we refer to the sumaclust manual.

    ${c_green}--vclust_params${c_reset}                   Additional parameters for vsearch cluster analysis. [default $params.vclust_params]
                                      For more information and options, we refer to the vsearch manual.

    ${c_green}--mmseqs_params${c_reset}                   Additional parameters for MMSeqs2 cluster analysis. [default $params.mmseqs_params]
                                      For more information and options, we refer to the MMSeqs2 manual.

    ${c_yellow}Computing options:${c_reset}
    ${c_green}--cores INT${c_reset}                       max cores per process for local use [default $params.cores]
    ${c_green}--max_cores INT${c_reset}                   max cores used on the machine for local use [default $params.max_cores]
    ${c_green}--memory INT${c_reset}                      max memory in GB for local use [default $params.memory]
    ${c_green}--output PATH${c_reset}                     name of the result folder [default $params.output]
    ${c_green}--permanentCacheDir PATH${c_reset}          location for auto-download data like databases [default $params.permanentCacheDir]
    ${c_green}--condaCacheDir PATH${c_reset}              location for storing the conda environments [default $params.condaCacheDir]
    ${c_green}--workdir PATH${c_reset}                    working directory for all intermediate results [default $params.workdir]

    ${c_yellow}Nextflow options:${c_reset}
    ${c_green}-with-report rep.html${c_reset}             cpu / ram usage (may cause errors)
    ${c_green}-with-dag chart.html${c_reset}              generates a flowchart for the process tree
    ${c_green}-with-timeline time.html${c_reset}          timeline (may cause errors)
    ${c_reset}____________________________________________________________________________________________

    """.stripIndent()
}