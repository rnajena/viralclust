/************************************************************************
* UPDATE_NCBI_METAINFO
*
* Downloads all viral genbank entries from the NCBI and pickles them
* to a Python3 dictionary for future annotation of sequences.
************************************************************************/

process update_ncbi_metainfo {
  label 'update_ncbi'

  input:

  output:

  script:
  """
  wget -N -q -P ${projectDir}/data/ ftp://ftp.ncbi.nih.gov/genbank/gbvrl*.seq.gz
  gunzip  ${projectDir}/data/*.gz
  python3 ${projectDir}/bin/ncbi_information_dump.py  ${projectDir}/data
  rm ${projectDir}/data/*seq ${projectDir}/data/*idx
  """

}