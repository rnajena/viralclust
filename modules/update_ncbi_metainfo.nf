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
  wget -c -N -q -P ${cacheDir}/ ftp://ftp.ncbi.nih.gov/genbank/gbvrl*.seq.gz
  wget -c -N -q -P ${cacheDir}/ ftp://ftp.ncbi.nih.gov/refseq/release/viral/*gbff.gz
  
  for FILE in ${cacheDir}/*.gz; do
    if [ -f \${FILE%.gz} ]; then
      continue
    fi
    gunzip \$FILE
  done
  echo "Starting Python Script"
  python3 ${projectDir}/bin/ncbi_information_dump.py  ${cacheDir}
  rm ${cacheDir}/*seq ${cacheDir}/*gbff ${cacheDir}/*idx
  """

}