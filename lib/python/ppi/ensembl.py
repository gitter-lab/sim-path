"""
Ensembl.org recommends the use of biomart to map HGNC, ENSG, ENST, ENSP
Currently these are implemented in the R script dl_ensembl_map
TODO add biomart python API functions here
"""

def map_hgnc_to_ensps(fh):
  """
  Return map from HGNC symbol to ENSP

  Parameters
  ----------
  fh : file-like 
    database file from dl_ensembl_map: 4 column TSV of HGNC, ENSG, ENST, ENSP

  Returns
  -------
  rv : dict
    mapping of HGNC symbol to a set of ENSPs
  """
  rv = {}
  for line in fh:
    line = line.rstrip()
    words = line.split('\t')
    hgnc_symbol = words[0]
    ensp_id = None
    if(len(words) >= 4):
      ensp_id = words[3]
    if ensp_id is not None:
      if hgnc_symbol not in rv:
        rv[hgnc_symbol] = set()
      rv[hgnc_symbol].add(ensp_id)
  return rv

def map_ensp_to_hgnc(fh):
  rv = {}
  for line in fh:
    line = line.rstrip()
    words = line.split('\t')
    hgnc_symbol = words[0]
    ensp_id = None
    if(len(words) >= 4):
      ensp_id = words[3]
    if ensp_id is not None:
      if ensp_id in rv:
        sys.stderr.write("Encountered ENSP id {} more than once\n".format(ensp_id))
      rv[ensp_id] = hgnc_symbol
  return rv

def apply_map(hgnc_to_ensp_map, gene_list_fh):
  """
  Parameters
  ----------
  hgnc_to_ensp_map : dict
    mapping of HGNC symbol to set of ENSPs

  gene_list_fh : file-like
    newline delimited file of HGNC symbols

  Returns
  -------
  rv : list
    sorted list of mapped ENSP identifiers
  """
  rv = set()
  for line in gene_list_fh:
    hgnc_symbol = line.rstrip()
    ensp_ids = hgnc_to_ensp_map.get(hgnc_symbol)
    if ensp_ids is None:
      sys.stderr.write("Unrecognized HGNC symbol: {}\n".format(hgnc_symbol))
    else:
      rv = rv.union(ensp_ids)
  rv = sorted(rv)
  return rv
