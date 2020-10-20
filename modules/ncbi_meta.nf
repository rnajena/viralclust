/************************************************************************
* NCBI_META
*
* Retrieve NCBI Meta information like taxonomy and collection date
************************************************************************/

process get_ncbi_meta {
  label 'ncbi_meta'

  input:
    path(sequences)

  output:
    path "ncbiMETA.pkl", emit: pkl_ncbi

  script:
  """
    python3 ${projectDir}/bin/get_ncbi_information.py "${sequences}"
  """
}