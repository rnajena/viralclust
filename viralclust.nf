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

if (params.hdbscan_help) {
  exit 0, hdbscanHelp()
}

if (params.mmseqs_help) {
  exit 0, mmseqsHelp()
}

if (params.sumaclust_help) {
  exit 0, sumaclustHelp()
}

if (params.cdhit_help) {
  exit 0, cdhitHelp()
}

if (params.vclust_help) {
  exit 0, vclustHelp()
}

if (params.hdbscan_off && params.cdhit_off && params.mmseqs_off && params.sumaclust_off && params.vclust_off) {
    exit 1, "ERROR: You switched off ALL cluster tools... Are you sure that makes sense?"
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
println "  $workflow.configFiles"
println "Path for database:"
println " $params.permanentCacheDir\u001B[0m"
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
    ${sw -> if (params.hdbscan_params != '') sw << "HDBSCAN parameters:     ${params.hdbscan_params}"}
    ${sw -> if (params.sumaclust_params != '') sw << "sumaclust parameters:     ${params.sumaclust_params}"}
    ${sw -> if (params.vclust_params != '') sw << "vclust parameters:     ${params.vclust_params}"}
    ${sw -> if (params.mmseqs_params != '') sw << "MMSeqs2 parameters:     ${params.mmseqs_params}"}


    """
    .stripIndent()

if (params.fasta) {

  sequences = Channel.fromPath(params.fasta)

  if (params.goi) {
    goi = Channel.fromPath(params.goi)
    include { sort_sequences as goi_sorter } from './modules/sortsequences'
    include { concat_goi } from './modules/remove_redundancy'
  }

  if (!params.sort_off) {
    include { sort_sequences } from './modules/sortsequences'
  }
  include { remove_redundancy } from './modules/remove_redundancy'

  include { cdhit } from './modules/cdhit'
  include { hdbscan } from './modules/hdbscan'
  include { sumaclust } from './modules/sumaclust'
  include { vclust } from './modules/vclust'
  include { mmseqs } from './modules/mmseqs'
  
  if (!params.sort_off) {
    include { reverseComp } from './modules/reverseComp'
  }
  include { mafft } from './modules/mafft'
  include { fasttree } from './modules/fasttree'
  include { nwdisplay } from './modules/nwutils'

  include { evaluate_cluster; merge_evaluation } from './modules/evaluate'
  include { generate_html } from './modules/generate_html'

  if (params.ncbi) {
    include { get_ncbi_meta } from './modules/ncbi_meta'
    ncbi_metainfo_ch = file("${params.permanentCacheDir}/ncbi_metainfo.pkl")
    if (! ncbi_metainfo_ch.exists() & ! params.update_ncbi) {
      include { update_ncbi_metainfo } from './modules/update_ncbi_metainfo'

      log.warn """ \
      No NCBI meta information database found.
      Database will be downloaded and stored for this and future runs!
      """.stripIndent()
    } else {
      p = "$params.permanentCacheDir/ncbi_metainfo.pkl"
      ncbiMeta = Channel.fromPath ( p )
    }
  }
}


if (params.update_ncbi) {
  include {update_ncbi_metainfo} from './modules/update_ncbi_metainfo'
}

workflow update_metadata {
  take:
    cacheDir
  main:
    update_ncbi_metainfo(cacheDir)
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
    if (!params.sort_off) {
        sort_sequences(sequences)
      if (params.goi) {
        goi_sorter(goi)
        goiSorted = goi_sorter.out.sort_result
      } else {
        goiSorted = 'NO FILE'
      }
      remove_redundancy(sort_sequences.out.sort_result)
    } else {
      remove_redundancy(sequences)
      goiSorted = 'NO FILE'
    }
    
    non_redundant_ch = remove_redundancy.out.nr_result
    if (params.goi) {
      concat_goi(remove_redundancy.out.nr_result, goiSorted)
      non_redundant_ch = concat_goi.out.nr_result
    }
  emit:
    non_redundant_ch
    goiSorted
}


workflow hdbscan_wf{
  take:
    fasta
    hdbscan_params
    goiSorted

  main:
    if (!params.hdbscan_off) {
      hdbscan(fasta, hdbscan_params, goiSorted)
      hdbscan_results = hdbscan.out[0]
    } else {
      hdbscan_results = Channel.from(['off', 'off', 'off'])
    }
  emit:
    hdbscan_results
}

workflow cd_hit_est_wf {
  take:
    fasta
    cdhit_params
    goiSorted
  
  main:
    if (!params.cdhit_off) {
      cdhit(fasta, cdhit_params, goiSorted)
    
      cdhit_results = cdhit.out[0]
    } else {
      cdhit_results = Channel.from(['off', 'off', 'off'])
    }

  emit:
    cdhit_results
}

workflow sumaclust_wf {
  take:
    fasta
    sumaclust_params
    goiSorted

  main: 
    if (!params.sumaclust_off) {
      sumaclust(fasta, sumaclust_params, goiSorted)
      sumaclust_results = sumaclust.out[0]
      } else {
        sumaclust_results = Channel.from(['off', 'off', 'off'])
      }
  
  emit:
    sumaclust_results

}

workflow vclust_wf {
  take:
    fasta
    vclust_params
    goiSorted

  main:
    if (!params.vclust_off) {
      vclust(fasta, vclust_params, goiSorted)
      vclust_results = vclust.out[0]
    } else {
      vclust_results = Channel.from(['off', 'off', 'off'])
    }

  emit:
    vclust_results
}

workflow mmseqs_wf {
  take:
    fasta
    mmseqs_params
    goiSorted

  main:
    if (!params.mmseqs_off) {
      mmseqs(fasta, mmseqs_params, goiSorted)
      mmseqs_results = mmseqs.out[0]
    } else {
      mmseqs_results = Channel.from(['off', 'off', 'off'])
    }
  
  emit:
    mmseqs_results
}

workflow clustering {

  take:
    non_redundant_ch
    goiSorted

  main:

    hdbscan_wf(non_redundant_ch, params.hdbscan_params, goiSorted)
    cd_hit_est_wf(non_redundant_ch, params.cdhit_params, goiSorted)
    sumaclust_wf(non_redundant_ch, params.sumaclust_params, goiSorted)
    vclust_wf(non_redundant_ch, params.vclust_params, goiSorted)
    mmseqs_wf(non_redundant_ch, params.mmseqs_params, goiSorted)

    results_channel = Channel.value('HDBSCAN').combine(hdbscan_wf.out)
                      .concat(Channel.value('cd-hit-est').combine(cd_hit_est_wf.out))
                      .concat(Channel.value('sumaclust').combine(sumaclust_wf.out))
                      .concat(Channel.value('MMseqs2').combine(mmseqs_wf.out))
                      .concat(Channel.value('vclust').combine(vclust_wf.out))
                      .filter { it[1] != 'off'}   
    if (!params.sort_off) {
      reverseComp(results_channel)
    }

  emit:
    results_channel
}

workflow phylo_wf {
  take:
    cluster_result

  main:
    if (params.eval || implicitEval) {
      mafft(cluster_result)
      fasttree(mafft.out.mafft_result)
      nwdisplay(fasttree.out.fasttree_result)

      tree_results = fasttree.out.fasttree_result
    } else {
      tree_results = cluster_result.flatten().filter(String).filter(~/[^.]*/).combine(Channel.from(['off']))
    }

  emit:
    tree_results  
}

workflow ncbi_wf {

  main:
    if (params.ncbi) {
      if (! ncbi_metainfo_ch.exists() & ! params.update_ncbi) {
        update_metadata("$params.permanentCacheDir")
        p = "$params.permanentCacheDir/ncbi_metainfo.pkl"
        ncbiMeta = Channel.fromPath ( p )
      }
      ncbiEval = Channel.from(eval_params).combine(ncbiMeta)
    } else {
      ncbiEval = Channel.from(['off']).combine(Channel.from(['off']))
    }

  emit:
    ncbiEval
}

workflow evaluation {
  take:
    cluster_result
    non_redundant_ch

  main:
    //test = cluster_result.flatten().filter(String).filter(~/[^.]*/)
    phylo_wf(cluster_result)
    ncbi_wf()
    
    eval_channel = cluster_result.combine(non_redundant_ch).
                    join(phylo_wf.out).
                    combine(ncbi_wf.out).view()

    evaluate_cluster(eval_channel)
    merge_evaluation(evaluate_cluster.out.eval_result.collect(), sequences)

    if (params.eval) {
    html_channel = cluster_result.
                    map{ it -> it[3] }.
                    collect()
    generate_html(html_channel)
    }


    // evaluate_cluster.out.warning.first().subscribe{

    //   log.warn """\

    //   ##########################################################
    //   NCBI meta information is older than 90 days.
    //   Please consider updating using the following command:
    //     nextflow run viralclust.nf --update_ncbi
    //   ##########################################################
    //   """.stripIndent()
    // }
}

workflow {

  if (params.update_ncbi) {
    update_metadata("$params.permanentCacheDir")
  }

  if (params.fasta) {
    preprocessing()
    clustering(preprocessing.out.non_redundant_ch, preprocessing.out.goiSorted)
    evaluation(clustering.out.results_channel, preprocessing.out.non_redundant_ch)
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

    ${c_green}--update_ncbi${c_reset}                     Downloads all current GenBank entries from the NCBI FTP server and processes the data to
                                      the databank stored at ${params.permanentCacheDir}.

    ${c_yellow}Cluster options:${c_reset}
    ${c_green}--cdhit_params${c_reset}                    Additional parameters for CD-HIT-EST cluster analysis. [default $params.cdhit_params]
                                      You can use ${c_green}nextflow run viralclust.nf --cdhit_help${c_reset}
                                      For more information and options, we refer to the CD-HIT manual.

    ${c_green}--hdbscan_params${c_reset}                  Additional parameters for HDBscan cluster analysis. [default $params.hdbscan_params]
                                      For more information and options, please use
                                      ${c_green}nextflow run viralclust.nf --hdbscan_help${c_reset}.

    ${c_green}--sumaclust_params${c_reset}                Additional parameters for sumaclust cluster analysis. [default $params.sumaclust_params]
                                      You can use ${c_green}nextflow run viralclust.nf --sumaclust_help${c_reset}.
                                      For more information and options, we refer to the sumaclust manual.

    ${c_green}--vclust_params${c_reset}                   Additional parameters for vsearch cluster analysis. [default $params.vclust_params]
                                      You can use ${c_green}nextflow run viralclust.nf --vclust_help${c_reset}
                                      For more information and options, we refer to the vsearch manual.

    ${c_green}--mmseqs_params${c_reset}                   Additional parameters for MMSeqs2 cluster analysis. [default $params.mmseqs_params]
                                      You can use ${c_green}nextflow run viralclust.nf --mmseqs_help${c_reset}
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

def hdbscanHelp() {
  c_reset = "\033[0m";
  c_red = "\033[1;31m"
  c_green = "\033[1;32m";
  c_yellow = "\033[0;33m";
  c_blue = "\033[0;34m";
  c_dim = "\033[2m";
  log.info """ \

  ____________________________________________________________________________________________

  This python program is part of ViralClust and takes several genome sequences
  from different viruses as an input.
  These will be clustered these sequences into groups (clades) based
  on their sequence similarity. For each clade, the centroid sequence is
  determined as representative genome, i.e. the sequence with the lowest
  distance to all other sequences of this clade.
____________________________________________________________________________________________

  Python Dependencies:
    docopt
    BioPython
    colorlog
    numpy
    scipy
    umap-learn
    hdbscan

  Contact:
    kevin.lamkiewicz@uni-jena.de


  Usage:
    hdbscan_virus.py [options] <inputSequences> [<genomeOfInterest>]

  Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
  --version                               Prints the version of viralClust and exits.
  -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]
  --metric METRIC                         Distance metric applied by UMAP (if applied) and HDBSCAN.
                                          The following are supported:
                                          'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                          'canberra', 'braycurtis',
                                          'mahalanobis', 'wminkowski', 'seuclidean',
                                          'cosine'.
                                          If an invalid metric is set, ViralClust will default back to 
                                          the cosine distance.
                                          [Default: cosine]

  --pca                                  Flag that determines whether (instead of UMAP) a PCA is used for dimension reduction. [Default: False]

  --neighbors NEIGHBORS                   Number of neighbors considered by UMAP to reduce the dimension space.
                                          Low numbers here mean focus on local structures within the data, whereas 
                                          larger numbers may loose fine details. [default: 50]
  --dThreshold dTHRESHOLD                 Sets the threshold for the minimum distance of two points in the low-dimensional space.
                                          Smaller thresholds are recommended for clustering and identifying finer topological structures
                                          in the data. [Default: 0.25]
  --dimension DIMENSION                   UMAP tries to find an embedding for the input data that can be represented by a low-dimensional space.
                                          This parameter tells UMAP how many dimensions should be used for the embedding. Lower numbers may result 
                                          in loss of information, whereas larger numbers will increase the runtime. [Default: 20]

  --clusterSize CLUSTERSIZE               This parameter forces HDBSCAN to form cluster with a size larger-equal to CLUSTERSIZE.
                                          Be aware that some data points (i.e. genomes) or even whole subcluster are considered as noise, if this parameter is set too high.
                                          E.g., if a very distinct viral genus has 40 genomes and the parameter is set to anything >40, HDBSCAN will not form
                                          the genus specific cluster. [Default: 5]
  --minSample MINSAMPLE                   Intuitively, this parameter declares how conservative clustering is performed. Higher values will lead 
                                          to more points considered noise, whereas a low value causes "edge-cases" to be grouped into a cluster.
                                          The default parameter is the same as CLUSTERSIZE. [Default: CLUSTERSIZE]
  ____________________________________________________________________________________________
  """.stripIndent()
}


def sumaclustHelp() {
  log.info """ \
------------------------------------------------------------
  SUMACLUST Version 1.0.31
------------------------------------------------------------
  Synopsis : star clustering of sequences.
  Usage: sumaclust [options] <dataset>
------------------------------------------------------------
  Options:
  -h       : [H]elp - print <this> help

  -l       : Reference sequence length is the shortest.

  -L       : Reference sequence length is the largest.

  -a       : Reference sequence length is the alignment length (default).

  -n       : Score is normalized by reference sequence length (default).

  -r       : Raw score, not normalized.

  -d       : Score is expressed in distance (default : score is expressed in similarity).

  -t ##.## : Score threshold for clustering. If the score is normalized and expressed in similarity (default),
            it is an identity, e.g. 0.95 for an identity of 95%. If the score is normalized
            and expressed in distance, it is (1.0 - identity), e.g. 0.05 for an identity of 95%.
            If the score is not normalized and expressed in similarity, it is the length of the
            Longest Common Subsequence. If the score is not normalized and expressed in distance,
            it is (reference length - LCS length).
            Only sequences with a similarity above ##.## with the center sequence of a cluster
            are assigned to that cluster. Default: 0.97.

  -e       : Exact option : A sequence is assigned to the cluster with the center sequence presenting the
            highest similarity score > threshold, as opposed to the default 'fast' option where a sequence is
            assigned to the first cluster found with a center sequence presenting a score > threshold.

  -R ##    : Maximum ratio between the counts of two sequences so that the less abundant one can be considered
            as a variant of the more abundant one. Default: 1.0.

  -p ##    : Multithreading with ## threads using openMP.

  -s ####  : Sorting by ####. Must be 'None' for no sorting, or a key in the fasta header of each sequence,
            except for the count that can be computed (default : sorting by count).

  -o       : Sorting is in ascending order (default : descending).

  -g       : n's are replaced with a's (default: sequences with n's are discarded).

  -B ###   : Output of the OTU table in BIOM format is activated, and written to file ###.

  -O ###   : Output of the OTU map (observation map) is activated, and written to file ###.

  -F ###   : Output in FASTA format is written to file ### instead of standard output.

  -f       : Output in FASTA format is deactivated.

------------------------------------------------------------
  Argument : the nucleotide dataset to cluster (or nothing
            if the standard input should be used).
------------------------------------------------------------
  http://metabarcoding.org/sumaclust
------------------------------------------------------------

  """.stripIndent()
}


def mmseqsHelp() {
  log.info """ \

  Clusters sequences by similarity in linear time. It groups similar sequences together based on user-specified criteria (max. E-value, seq. id., min. coverage,...).
  By Martin Steinegger <martin.steinegger@mpibpc.mpg.de>

Options:
  Prefilter:
    --comp-bias-corr INT            correct for locally biased amino acid composition (range 0-1) [1]
    --add-self-matches              artificially add entries of queries with themselves (for clustering)
    --alph-size INT                 alphabet size (range 2-21) [21]
    --mask INT                      mask sequences in k-mer stage 0: w/o low complexity masking, 1: with low complexity masking [1]
    --mask-lower-case INT           lowercase letters will be excluded from k-mer search 0: include region, 1: exclude region [0]
    -k INT                          k-mer size in the range (0: set automatically to optimum) [0]
    --split-memory-limit BYTE       Set max memory per split. E.g. 800B, 5K, 10M, 1G. Defaults (0) to all available system memory. [0]

  Align:
    -a                              add backtrace string (convert to alignments with mmseqs convertalis utility)
    --alignment-mode INT            How to compute the alignment: 0: automatic; 1: only score and end_pos; 2: also start_pos and cov; 3: also seq.id; 4: only ungapped alignment [0]
    -e FLOAT                        list matches below this E-value (range 0.0-inf) [0.001]
    --min-seq-id FLOAT              list matches above this sequence identity (for clustering) (range 0.0-1.0) [0.000]
    --min-aln-len INT               minimum alignment length (range 0-INT_MAX) [0]
    --seq-id-mode INT               0: alignment length 1: shorter, 2: longer sequence [0]
    --alt-ali INT                   Show up to this many alternative alignments [0]
    -c FLOAT                        list matches above this fraction of aligned (covered) residues (see --cov-mode) [0.800]
    --cov-mode INT                  0: coverage of query and target, 1: coverage of target, 2: coverage of query 3: target seq. length needs to be at least x% of query length, 4: query seq. length needs to be at least x% of target length 5: short seq. needs to be at least x% of the other seq. length [0]
    --realign                       compute more conservative, shorter alignments (scores and E-values not changed)
    --max-rejected INT              maximum rejected alignments before alignment calculation for a query is aborted [2147483647]
    --max-accept INT                maximum accepted alignments before alignment calculation for a query is stopped [2147483647]
    --score-bias FLOAT              Score bias when computing the SW alignment (in bits) [0.000]
    --gap-open INT                  Gap open cost [11]
    --gap-extend INT                Gap extension cost [1]

  Clust:
    --cluster-mode INT              0: Setcover, 1: connected component, 2: Greedy clustering by sequence length  3: Greedy clustering by sequence length (low mem) [0]
    --max-iterations INT            maximum depth of breadth first search in connected component [1000]
    --similarity-type INT           type of score used for clustering (range 1,2). 1=alignment score. 2=sequence identity  [2]

  Kmermatcher:
    --kmer-per-seq INT              kmer per sequence [21]
    --adjust-kmer-len               adjust k-mer length based on specificity (only for nucleotides)
    --hash-shift INT                Shift k-mer hash [5]
    --include-only-extendable       Include only extendable
    --skip-n-repeat-kmer INT        Skip sequence with >= n exact repeating k-mers [0]

  Profile:
    --pca FLOAT                     pseudo count admixture strength [1.000]
    --pcb FLOAT                     pseudo counts: Neff at half of maximum admixture (range 0.0-inf) [1.500]

  Misc:
    --rescore-mode INT              Rescore diagonal with: 0: Hamming distance, 1: local alignment (score only), 2: local alignment, 3: global alignment or 4: longest alignment fullfilling window quality criterion [0]
    --remove-tmp-files 0            Delete temporary files [1, set to 0 to disable]
    --dont-split-seq-by-len 0       Dont split sequences by --max-seq-len [1, set to 0 to disable]
    --dbtype INT                    Database type 0: auto, 1: amino acid 2: nucleotides [0]
    --dont-shuffle 0                Do not shuffle input database [1, set to 0 to disable]
    --id-offset INT                 numeric ids in index file are offset by this value  [0]

  Common:
    --threads INT                   number of cores used for the computation (uses all cores by default) [8]
    --compressed INT                write results in compressed format [0]
    -v INT                          verbosity level: 0=nothing, 1: +errors, 2: +warnings, 3: +info [3]
    --sub-mat MAT                   amino acid substitution matrix file [nucl:nucleotide.out,aa:blosum62.out]
    --max-seq-len INT               maximum sequence length (range 1-32768]) [65535]
    --db-load-mode INT              Database preload mode 0: auto, 1: fread, 2: mmap, 3: mmap+touch [0]
    --force-reuse                   reuse tmp file in tmp/latest folder ignoring parameters and git version change
    --mpi-runner STR                Use MPI on compute grid with this MPI command (e.g. "mpirun -np 42") []

  Expert:
    --kmer-per-seq-scale FLOAT      scale kmer per sequence based on sequence length as kmer-per-seq val + scale x seqlen [0.000]
    --filter-hits                   filter hits by seq.id. and coverage
    --sort-results INT              Sort results: 0: no sorting, 1: sort by evalue (Alignment) or seq.id. (Hamming) [0]

  """.stripIndent()
}


def vclustHelp() {

  log.info """ \

  General options
  --bzip2_decompress          decompress input with bzip2 (required if pipe)
  --fasta_width INT           width of FASTA seq lines, 0 for no wrap (80)
  --gzip_decompress           decompress input with gzip (required if pipe)
  --help | -h                 display help information
  --log FILENAME              write messages, timing and memory info to file
  --maxseqlength INT          maximum sequence length (50000)
  --minseqlength INT          min seq length (clust/derep/search: 32, other:1)
  --no_progress               do not show progress indicator
  --notrunclabels             do not truncate labels at first space
  --quiet                     output just warnings and fatal errors to stderr
  --threads INT               number of threads to use, zero for all cores (0)
  --version | -v              display version information


  Clustering
  --cluster_fast FILENAME     cluster sequences after sorting by length
  --cluster_size FILENAME     cluster sequences after sorting by abundance
  --cluster_smallmem FILENAME cluster already sorted sequences (see -usersort)
  --cluster_unoise FILENAME   denoise Illumina amplicon reads
  Parameters (most searching options also apply)
  --cons_truncate             do not ignore terminal gaps in MSA for consensus
  --id REAL                   reject if identity lower, accepted values: 0-1.0
  --iddef INT                 id definition, 0-4=CD-HIT,all,int,MBL,BLAST (2)
  --qmask none|dust|soft      mask seqs with dust, soft or no method (dust)
  --sizein                    propagate abundance annotation from input
  --strand plus|both          cluster using plus or both strands (plus)
  --usersort                  indicate sequences not pre-sorted by length
  --minsize INT               minimum abundance (unoise only) (8)
  --unoise_alpha REAL         alpha parameter (unoise only) (2.0)
  Output
  --biomout FILENAME          filename for OTU table output in biom 1.0 format
  --centroids FILENAME        output centroid sequences to FASTA file
  --clusterout_id             add cluster id info to consout and profile files
  --clusterout_sort           order msaout, consout, profile by decr abundance
  --clusters STRING           output each cluster to a separate FASTA file
  --consout FILENAME          output cluster consensus sequences to FASTA file
  --mothur_shared_out FN      filename for OTU table output in mothur format
  --msaout FILENAME           output multiple seq. alignments to FASTA file
  --otutabout FILENAME        filename for OTU table output in classic format
  --profile FILENAME          output sequence profile of each cluster to file
  --relabel STRING            relabel centroids with this prefix string
  --relabel_keep              keep the old label after the new when relabelling
  --relabel_md5               relabel with md5 digest of normalized sequence
  --relabel_self              relabel with the sequence itself as label
  --relabel_sha1              relabel with sha1 digest of normalized sequence
  --sizeorder                 sort accepted centroids by abundance, AGC
  --sizeout                   write cluster abundances to centroid file
  --uc FILENAME               specify filename for UCLUST-like output
  --xsize                     strip abundance information in output

  """.stripIndent()
}

def cdhitHelp() {
  log.info """\
        ====== CD-HIT version 4.6 (built on Aug 29 2016) ======

Usage: cd-hit-est [Options]

Options

    -i	input filename in fasta format, required
    -o	output filename, required
    -c	sequence identity threshold, default 0.9
    this is the default cd-hit's "global sequence identity" calculated as:
    number of identical amino acids in alignment
    divided by the full length of the shorter sequence
    -G	use global sequence identity, default 1
    if set to 0, then use local sequence identity, calculated as :
    number of identical amino acids in alignment
    divided by the length of the alignment
    NOTE!!! don't use -G 0 unless you use alignment coverage controls
    see options -aL, -AL, -aS, -AS
    -b	band_width of alignment, default 20
    -M	memory limit (in MB) for the program, default 800; 0 for unlimitted;
    -T	number of threads, default 1; with 0, all CPUs will be used
    -n	word_length, default 10, see user's guide for choosing it
    -l	length of throw_away_sequences, default 10
    -d	length of description in .clstr file, default 20
    if set to 0, it takes the fasta defline and stops at first space
    -s	length difference cutoff, default 0.0
    if set to 0.9, the shorter sequences need to be
    at least 90% length of the representative of the cluster
    -S	length difference cutoff in amino acid, default 999999
    if set to 60, the length difference between the shorter sequences
    and the representative of the cluster can not be bigger than 60
    -aL	alignment coverage for the longer sequence, default 0.0
    if set to 0.9, the alignment must covers 90% of the sequence
    -AL	alignment coverage control for the longer sequence, default 99999999
    if set to 60, and the length of the sequence is 400,
    then the alignment must be >= 340 (400-60) residues
    -aS	alignment coverage for the shorter sequence, default 0.0
    if set to 0.9, the alignment must covers 90% of the sequence
    -AS	alignment coverage control for the shorter sequence, default 99999999
    if set to 60, and the length of the sequence is 400,
    then the alignment must be >= 340 (400-60) residues
    -A	minimal alignment coverage control for the both sequences, default 0
    alignment must cover >= this value for both sequences
    -uL	maximum unmatched percentage for the longer sequence, default 1.0
    if set to 0.1, the unmatched region (excluding leading and tailing gaps)
    must not be more than 10% of the sequence
    -uS	maximum unmatched percentage for the shorter sequence, default 1.0
    if set to 0.1, the unmatched region (excluding leading and tailing gaps)
    must not be more than 10% of the sequence
    -U	maximum unmatched length, default 99999999
    if set to 10, the unmatched region (excluding leading and tailing gaps)
    must not be more than 10 bases
    -B	1 or 0, default 0, by default, sequences are stored in RAM
    if set to 1, sequence are stored on hard drive
    it is recommended to use -B 1 for huge databases
    -p	1 or 0, default 0
    if set to 1, print alignment overlap in .clstr file
    -g	1 or 0, default 0
    by cd-hit's default algorithm, a sequence is clustered to the first
    cluster that meet the threshold (fast cluster). If set to 1, the program
    will cluster it into the most similar cluster that meet the threshold
    (accurate but slow mode)
    but either 1 or 0 won't change the representatives of final clusters
    -r	1 or 0, default 1, by default do both +/+ & +/- alignments
    if set to 0, only +/+ strand alignment
    -mask	masking letters (e.g. -mask NX, to mask out both 'N' and 'X')
    -match	matching score, default 2 (1 for T-U and N-N)
    -mismatch	mismatching score, default -2
    -gap	gap opening score, default -6
    -gap-ext	gap extension score, default -1
    -bak	write backup cluster file (1 or 0, default 0)
    -h	print this help

    Questions, bugs, contact Limin Fu at l2fu@ucsd.edu, or Weizhong Li at liwz@sdsc.edu
    For updated versions and information, please visit: http://cd-hit.org

    cd-hit web server is also available from http://cd-hit.org

    If you find cd-hit useful, please kindly cite:

    "Clustering of highly homologous sequences to reduce thesize of large protein database", Weizhong Li, Lukasz Jaroszewski & Adam Godzik. Bioinformatics, (2001) 17:282-283
    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences", Weizhong Li & Adam Godzik. Bioinformatics, (2006) 22:1658-1659
  """.stripIndent()
}
