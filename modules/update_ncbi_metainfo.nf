/************************************************************************
* UPDATE_NCBI_METAINFO
*
* Downloads all viral genbank entries from the NCBI and pickles them
* to a Python3 dictionary for future annotation of sequences.
************************************************************************/

process update_ncbi_metainfo {
  label 'update_ncbi'

  input:
    path(cacheDir)

  output:

  script:
  """
  wget -N -q -P ${projectDir}/${cacheDir}/ ftp://ftp.ncbi.nih.gov/genbank/gbvrl*.seq.gz
  gunzip  ${projectDir}/${cacheDir}/*.gz
  python3 ${projectDir}/bin/ncbi_information_dump.py  ${projectDir}/${cacheDir}
  rm ${projectDir}/${cacheDir}/*seq ${projectDir}/${cacheDir}/*idx
  """

}