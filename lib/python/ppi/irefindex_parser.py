import os.path
import pickle
import re
import numpy as np
import copy
import scipy.sparse as sp
import sys
import itertools as it

# for write_lpr_hist only
import scipy.stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# TODO function to list allowed identifier types
# NOTE we are only interested in homo sapiens at this time
TAX_ID = "taxid:9606(Homo sapiens)"
TAX_A = "taxa"
TAX_B = "taxb"
ALIAS_A = "aliasA"
ALIAS_B = "aliasB"
CHECKSUM = "checksum"
CHECKSUM_A = "Checksum_A"
CHECKSUM_B = "Checksum_B"
IROGID_A = "irogida"
IROGID_B = "irogidb"
HGNCS_KEY = "hgncs"
CONFIDENCE = "confidence"
LPR = "lpr" # "lowest publication record"?
KNOWN_ALIASES = ["crogid", "entrezgene/locuslink", "genbank_protein_gi", "hgnc", "icrogid", "irogid", "pdb", "refseq", "rogid", "uniprotkb"] # see get_alias_names

def write_abc(ifh, ofh, id_type='rogid'):
  """
  Translate irefindex database file to a graph file in "ABC" format, a tab-delimited graph file format
    <interactor A>	<interactor B>	<weight>

  Parameters
  ----------
  ifh : io-like
    input file handle containing irefindex data
  
  ofh : io-like
    output file handle to write ABC file to

  Returns
  -------
  num_edges : int
    Number of interactions written
  """
  num_edges = 0

  try:
    ifh.seek(0)
  except AttributeError:
    raise ArgumentError("ifh must implement seek")

  # first pass through ifh to determine range of lpr values
  header = parse_header(ifh)
  ifh.seek(0)
  min_lpr, max_lpr = get_lpr_extrema(ifh)
  for line in ifh:
    line = line.rstrip()
    datum = line.split("\t")
    lpr = get_lpr(header, datum)
    if lpr is not None:
      if(lpr < min_lpr):
        min_lpr = lpr
      if(lpr > max_lpr):
        max_lpr = lpr

  def write_line(id_a, id_b, weight):
    abc = [id_a, id_b, str(weight)]
    out_line = "\t".join(abc) + "\n"
    ofh.write(out_line)

  ifh.seek(0) # if havent raised previously, can do
  parse_header(ifh) # advance past header
  for line in ifh:
    line = line.rstrip()
    datum = line.split("\t")
    lpr = get_lpr(header, datum)
    if lpr is None:
      # then assign worst confidence score
      lpr = max_lpr
    weight = lpr_to_conf(min_lpr, max_lpr, lpr)
    if(id_type == 'rogid'):
      id_a, id_b = get_checksums(header, datum)
      write_line(id_a, id_b, weight)
      num_edges += 1
    elif(id_type == 'hgnc'):
      ids_a, ids_b = get_hgnc_pair(header, datum)
      missing_hgnc = (len(ids_a) == 0 or len(ids_b) == 0)
      if(not missing_hgnc):
        for id_a, id_b in it.product(ids_a, ids_b):
          write_line(id_a, id_b, weight)
          num_edges += 1
    else:
      raise Exception('unrecognized id type')

  return num_edges

def lpr_to_conf(min_lpr, max_lpr, lpr):
  # linearly transform lpr to a confidence score in 0 to 1; where min_lpr is 1 and max_lpr is 0
  # TODO other transformation choices?
  length = max_lpr - min_lpr
  conf = float(lpr - min_lpr) / length
  return conf

def write_lpr_hist(ifh, ofp):
  """
  Parse irefindex and create a histogram of lpr values
  """
  header = parse_header(ifh)
  lprs = []
  n_none = 0
  for line in ifh:
    line = line.rstrip()
    datum = line.split("\t")
    lpr = get_lpr(header, datum)
    if lpr is None:
      n_none += 1
    else:
      lprs.append(lpr)
  lprs = np.array(lprs)

  n_bins = 50
  n, bins, patches = plt.hist(lprs, bins=n_bins)

  plt.title("LPR values; n_none: {}".format(n_none))
  plt.xlabel("LPR")
  plt.yscale("log")
  plt.ylabel("Frequency")
  plt.savefig(ofp)

def get_lpr(header, datum):
  """
  Get LPR as an int or return None
  """
  conf = datum[header[CONFIDENCE]]
  avp_strs = conf.split("|") # attribute-value pairs
  lpr = None
  for avp_str in avp_strs:
    avp = avp_str.split(":")
    if(avp[0] == LPR):
      lpr = int(avp[1])
  return lpr

def get_lpr_extrema(fh):
  header = parse_header(fh)
  min_lpr = np.inf
  max_lpr = 0
  for line in fh:
    line = line.rstrip()
    datum = line.split("\t")
    lpr = get_lpr(header, datum)
    if lpr is not None:
      if(lpr < min_lpr):
        min_lpr = lpr
      if(lpr > max_lpr):
        max_lpr = lpr
  return min_lpr, max_lpr

# TODO deprecated
def get_row_vec(ids, alias_map, aliases):
  """Analogous to a somatic mutation profile in Hofree: Construct a column vector of the same 
  size as a PPI network where each of the interactors in @aliases@ has a 1 in its associated
  entry

  Parameters
  ----------
  ids : list
    Given a PPI network represented as an adjacency matrix where the (irefindex uid) identifiers for the interactor in each row are provided in a list of the same size, @ids@, 

  alias_map : dict
    See rekey_aliases

  aliases: list
    aliases of interest

  Returns
  -------
  row_vec : shape (1,n) scipy sparse array

  missing : list
    subset of aliases that cannot be found
  """
  row_vec = np.zeros(len(ids))
  missing = []

  # TODO speed up with via better search of ids
  for alias in aliases:
    if alias in alias_map:
      checksum = alias_map[alias][CHECKSUM]
      index = ids.index(checksum)
      row_vec[index] = 1
    else:
      missing.append(alias)

  row_vec_sparse = sp.csc_matrix(row_vec)
  return (row_vec_sparse, missing)

def get_row_vec_alt(all_ids, ids):
  row_vec = np.zeros(len(all_ids))
  missing = []
  for id in ids:
    try:
      index = all_ids.index(id)
      row_vec[index] = 1
    except ValueError:
      missing.append(id)

  row_vec_sparse = sp.csc_matrix(row_vec)
  return (row_vec_sparse, missing)

def parse_header(ifh):
  """Modify ifh by advancing past the first line and parsing the header
  @return dict mapping string of field name to its index
  """
  field_map = {}
  i = 0
  # ignore leading "#"
  for field in ifh.readline()[1:].split("\t"):
    field_map[field] = i
    i += 1
  return field_map

def get_hgnc(header, data):
  """Using header, extract hgnc names from interaction data
  @param header @see parse_header
  @param data a data row from iRefIndex that has been split on tabs
  @return list of all HGNC identifiers for this interaction
  """
  rv = []
  hgnc_regexp = re.compile('^hgnc:(.*)')
  alias_as = data[header[ALIAS_A]].split("|")
  for alias_a in alias_as:
    match_data = hgnc_regexp.match(alias_a)
    if(match_data is not None):
      rv.append(match_data.group(1).upper())

  alias_bs = data[header[ALIAS_B]].split("|")
  for alias_b in alias_bs:
    match_data = hgnc_regexp.match(alias_b)
    if(match_data is not None):
      rv.append(match_data.group(1).upper())

  return rv

def get_hgnc_pair(header, datum):
  """
  See get_hgnc. While get_hgnc gets all hgnc identifiers, this function
  keeps the HGNC identifiers for interactor A separate from those for
  interactor B.
  """
  return get_id_pair(header, datum, id_type='hgnc')

def get_id_pair(header, datum, id_type='hgnc'):
  """
  Return a list of identifiers of type <id_type> for each interactor A and B in the current record <datum>

  Parameters
  ----------
  id_type : string
    member of KNOWN_ALIASES
  """
  all_ids_a, all_ids_b = get_identifiers(header, datum)
  ids_a = all_ids_a.get(id_type, [])
  ids_b = all_ids_b.get(id_type, [])
  rv = (list(ids_a), list(ids_b))
  return rv

def get_checksums(header, datum):
  """See get_identifiers

  Returns
  -------
  checksum_a : 2-tuple
    unique identifier for interactor A

  checksum_b : 2-tuple
    see checksum_a
  """
  a_dict = _get_identifiers(header, datum, [CHECKSUM_A])
  b_dict = _get_identifiers(header, datum, [CHECKSUM_B])
  if len(a_dict) > 1 or len(b_dict) > 1:
    raise RuntimeError("did not expect multiple checksum TYPES for datum {}".format(datum))

  checksum_a = checksum_b = None
  if(sys.version_info[0] == 3):
    # TODO modify to match else
    raise RuntimeError("not implemented")
    checksum_a = list(a_dict.items())[0]
    checksum_b = list(b_dict.items())[0]
  else:
    list_a = list(a_dict.items()[0][1])
    list_b = list(b_dict.items()[0][1])
    if len(list_a) > 1 or len(list_b) > 1:
      print(list_a)
      print(list_b)
      raise RuntimeError("did not expect multiple checksum VALUES for datum {}".format(datum))
    checksum_a = list_a[0]
    checksum_b = list_b[0]

  return (checksum_a, checksum_b)

def get_identifiers(header, data):
  """See get_hgnc, extract all identifiers not just HGNC

  Parameters
  ----------
  header : dict
    map of fields from irefindex header line to their index

  data : list
    a data row in irefindex with each field as an element in the list

  Returns
  -------
  ids_a : dict
    map of identifier name to a set of values for that identifier name

  ids_b : dict
    see ids_a

  TODO
  ----
  there may be multiple identifiers from the same database as in multiple uniprotkb identifiers:
    'uniprotkb': 'IFIT2_HUMAN', 'uniprotkb' 'P09913'
  """
  # TODO update functions using this one to use the values as a set rather than a single string
  fields_a = ["uidA", "altA", "aliasA"]
  fields_b = ["uidB", "altB", "aliasB"]

  ids_a = _get_identifiers(header, data, fields_a)
  ids_b = _get_identifiers(header, data, fields_b)

  return (ids_a, ids_b)

def _get_identifiers(header, data, fields):
  """See get_identifiers

  Returns
  -------
  rv : dict
    map from field in <fields> to a set of identifiers
  """
  rv = {}
  for field in fields:
    pipe_str = data[header[field]]
    id_value_pairs = pipe_str.split("|")
    for id_value_pair in id_value_pairs:
      # some identifiers have ":" in the value; but other dont like irogid
      try:
        i = id_value_pair.index(":")
        id = id_value_pair[0:i]
        value = id_value_pair[i+1:]
        if(id not in rv):
          rv[id] = set()
        rv[id].add(value)
      except ValueError:
        # then no ":" in the value and id_value_pair is just a value
        if(field not in rv):
          rv[field] = set()
        rv[field].add(id_value_pair)
  return rv

def get_hgnc_universe(ifh):
  """Extract all HGNC identifiers from iRefIndex data file @ifh@ and return them as a set"""
  header = parse_header(ifh)
  hgnc_set = set()
  for line in ifh:
    data = line.split("\t")
    hgncs = get_hgnc(header, data)
    for hgnc in hgncs:
      hgnc_set.add(hgnc)
  return hgnc_set

def save_hgnc(ifh, ofh):
  """Extract all HGNC identifiers from iRefIndex data file @ifh@ and save them in @ofh@
  """
  hgnc_set = get_hgnc_universe(ifh)
  hgnc_list = list(hgnc_set)
  hgnc_list = sorted(hgnc_list)
  for hgnc in hgnc_list:
    ofh.write(hgnc + "\n")

def verify_hgnc(all_hgncs, cand_hgncs):
  """Verify that the HGNCSs in <cand_hgncs> are found in <all_hgncs>. Partition
  <cand_hgncs> into two sets: those that are found in all_hgncs and those that are not.

  Parameters
  ----------
  all_hgncs : list
    universe of HGNCs we are considering

  cand_hgncs : list
    list of HGNCs (from another source) that should be verified

  Returns
  -------
  (found, missing) : (list, list)
  """
  all_hgncs = set(all_hgncs)
  found = filter(lambda x: x in all_hgncs, cand_hgncs)
  missing = filter(lambda x: x not in all_hgncs, cand_hgncs)
  return (found, missing)

def is_human_human(header, datum):
  """Return true if interaction data is a human-human interaction
  """
  return (datum[header[TAX_A]] == TAX_ID and datum[header[TAX_B]] == TAX_ID)

def map_rogid_to_hgnc(fh):
  """
  Parameters
  ----------
  fh : io-like
    irefindex file handle

  Returns
  -------
  rv : dict
    mapping of rogid to a set of hgnc identifiers
  """
  rv = {}
  rogid_to_aliases = get_aliases(fh)
  for rogid, alias_to_id in rogid_to_aliases.items():
    hgnc = None
    if 'hgnc' in alias_to_id:
      hgnc = alias_to_id['hgnc']
    if rogid in rv and rv[rogid] != hgnc:
      raise Exception("rogid {} has at least 2 hgnc identifiers {} and {}".format(rogid, hgnc, rv[rogid]))
    rv[rogid] = hgnc
  return rv

def map_hgnc_to_rogids(rogid_to_hgnc):
  rv = {}
  for rogid, hgnc_set in rogid_to_hgnc.items():
    if hgnc_set is not None:
      for hgnc in hgnc_set:
        if(hgnc not in rv):
          rv[hgnc] = set()
        rv[hgnc].add(rogid)
  return rv

def get_aliases(fh):
  """Map irefindex unique identifiers (checksum of taxon and amino acid sequence, not "uidA")
  to a dictionary of their aliases

  Parameters
  -----------
  fh : file handle
    irefindex data file

  Returns
  -------
  alias_map : dict
    map of uid to aliases; aliases a dictionary mapping identifier name to a set of identifier values;
    see get_identifiers
  """
  try:
    fh.seek(0)
  except AttributeError:
    pass

  alias_map = {}
  header = parse_header(fh)
  for line in fh:
    line = line.strip()
    datum = line.split("\t")

    # include only human-human interactions
    if (is_human_human(header, datum)):
      # get identifiers from interaction data
      checksum_a, checksum_b = get_checksums(header, datum)
      ids_a, ids_b = get_identifiers(header, datum)

      # add ids to alias_map, ensuring that all irefindex identifiers are consistent
      if checksum_a in alias_map:
        assert ids_a == alias_map[checksum_a], "inconsistent identifiers for checksum {}: {} ; {}".format(checksum_a, ids_a, alias_map[checksum_a])
      else:
        alias_map[checksum_a] = ids_a
      if checksum_b in alias_map:
        assert ids_b == alias_map[checksum_b], "inconsistent identifiers for checksum {}: {} ; {}".format(checksum_b, ids_b, alias_map[checksum_b])
      else:
        alias_map[checksum_b] = ids_b

  return alias_map

def get_alias_names(alias_map):
  """
  Return a union of all alias keys
  See get_aliases
  """
  rv = set()
  for id_v, aliases in alias_map.iteritems():
    rv = rv.union(aliases.keys())
  return list(rv)

def rekey_aliases(alias_map, rekey):
  """Transform alias_map into a new dictionary mapping the identifier type @rekey@ to its aliases

  Parameters
  ----------
  alias_map : dict
    See get_aliases

  rekey : hashable
    e.g. 'hgnc'

  Returns
  -------
  new_alias_map : dict
  """
  new_alias_map = {}
  for k, alias_to_values in alias_map.iteritems():
    alias_to_values = copy.deepcopy(alias_to_values)
    if rekey in alias_to_values:
      revalues = alias_to_values[rekey]
      # TODO what if key uid already exists
      alias_to_values[CHECKSUM] = k
      # TODO upper() only needed for keys like HGNC but not others
      for revalue in revalues:
        new_alias_map[revalue.upper()] = alias_to_values
  return new_alias_map

def get_alias_for_id(alias_map, id_v, id_type):
  """
  Return aliases_for_id, a dictionary mapping id_type to id_value or None if the id is not in irefindex
  """
  re_alias_map = rekey_aliases(alias_map, id_type)
  aliases_for_id = re_alias_map.get(id_v)
  return aliases_for_id

def translate_id(alias_map, type_src, type_dest, id_src):
  """Return <id_dest>, an identifier alias of type <type_dest> for the identifier given by <id_src>
  """
  id_dest = None
  aliases = get_alias_for_id(alias_map, id_src, type_src)
  if aliases is not None:
    id_dest = aliases.get(type_dest)
  return id_dest

def parse_irefindex(path, n_data=5):
  """Parse first @n_data@ data lines irefindex data file at @path@
  @return field_map, data where @field_map@ is a map of field name to index number
    and @data@ is an array of arrays with inner array the result of splitting a 
    data line on tab characters
  """
  field_map = {}
  data = []
  with open(path) as fh:
    header_line = fh.readline()
    i = 0
    for field in header_line[1:].split("\t"):
      field_map[field] = i
      i += 1

    i = 0
    for line in fh:
      data.append(line.split("\t"))
      i += 1
      if (i > n_data):
        break
  return field_map, data

def pretty_print(header, datum):
  """Print an irefindex data row for visual inspection
  """
  items = sorted(header.items(), lambda x,y: cmp(x[1], y[1]))
  for field, index in items:
    print("{}: {}".format(field, datum[header[field]]))

def map_id_to_index(ids_fh):
  """Returns a map of a rogid to its index in the diffusion array
  """
  try:
    ids_fh.seek(0)
  except AttributeError:
    pass

  ids_list = pickle.load(ids_fh)
  rv = {}
  i = 0
  for id_v in ids_list:
    rv[id_v] = i
    i += 1

  return rv

def rekey_id_map(alias_map, id_map, id_type):
  """Transform id_map into a new dictionary mapping identifier of type <id_type> to its index

  Parameters
  ----------
  alias_map : dict
    See get_aliases

  id_map : dict
    See map_id_to_index

  id_type : string
    See get_alias_names; id_type must be one of the values returned by that function

  Returns
  -------
  (re_id_map, missing_aliases) : (dict, list)
    missing_aliases is a list 
  """
  re_id_map = {}
  missing_aliases = []
  for k, index in id_map.iteritems():
    aliases = alias_map.get(k)
    if aliases is None:
      # require that id_map contains only known rogid values
      raise 'cannot find key {} in alias_map'.format(k)

    id_v = aliases.get(id_type)
    if id_v is None:
      # then there is no alias of id_type
      missing_aliases.append(k)
    else:
      # note: some rogid map to same hgnc, for example
      re_id_map[id_v] = index

  return (re_id_map, missing_aliases)

def filter_ids(alias_map, id_type, ids):
  """Partition <ids> into two lists. The first list contains the ids that are found in the
  network. The second list contains the ids that are not found in the network.

  Parameters
  ----------
  alias_map : dict
    See get_aliases

  id_type : string
    See get_alias_names; id_type must be one of the values returned by that function

  ids : list
    list of identifiers with type <id_type>

  Returns
  -------
  (found, missing) : (list, list)
  """
  re_alias_map = rekey_aliases(alias_map, id_type)
  # TODO what if id_type is bad
  found = []
  missing = []
  ids = map(lambda x: x.upper(), ids)
  for id in ids:
    if(alias_map.has_key(id)):
      found.append(id)
    else:
      missing.append(id)

  return (found, missing)

def lookup_loadings(loadings_arr, ids, alias_map, hgncs):
  """Return array parallel to hgncs where each element contains the loading factor 
  coefficients for that gene

  Parameters
  ----------
  loadings_arr : numpy array
    full loading factor matrix, "W" in X = W * H where X is (n_gene x n_list),
    W is (n_gene x n_factor)
    H is (n_factor x n_list)

  ids : list
    each element is a rogid identifier that is associated with index i in <loadings_arr>

  alias_map : dict
    See get_aliases; mapping of rogid to hgnc (and other aliases)

  hgncs : list
    list of hgnc identifiers to lookup loadings for

  Returns
  -------
  hgnc_loadings : list
    zipped list of hgnc identifier with numpy row vector of latent factor loadings associated with it
  """
  hgnc_id_map, missing_aliases = rekey_id_map(alias_map, ids, 'hgnc')
  indexes = map(lambda x: hgnc_id_map.get(x), hgncs)
  if(loadings_arr.shape[0] < loadings_arr.shape[1]):
    loadings_arr = loadings_arr.transpose()
  loadings = map(lambda x: loadings_arr[x, :] if x is not None else None, indexes)
  hgnc_loadings = zip(hgncs, loadings)
  return hgnc_loadings

def save_loadings(fh, hgnc_loadings):
  """Save hgnc_loadings to a file
  """
  # first, determine the number of latent factors from the first
  # hgnc that exists in irefindex
  n_factors = None
  for hgnc_loading in hgnc_loadings:
    hgnc, loadings = hgnc_loading
    if loadings is not None:
      n_factors = loadings.shape[0]
      break

  if n_factors is None:
    raise "unable to determine n_factors because hgnc_loadings contains no known hgnc identifiers"

  for hgnc_loading in hgnc_loadings:
    hgnc, loadings = hgnc_loading
    if loadings is None:
      # then, use found n_factors to normalize missing data
      loadings = np.zeros(n_factors)

    # next, save the loadings
    line = ",".join([hgnc] + map(str, loadings))
    fh.write(line + "\n")

def tmp(loadings, hgnc_loadings):
  """Determine which hgnc loadings exceed a given percentile of loading scores
  """
  # NOTE yields only those which are in the original gene list
  # TODO
  vs = map(lambda x: x[1], hgnc_loadings)
  vs2 = map(lambda x: x if x is not None else np.array((0,0,0,0,0)), vs)
  vs3 = map(lambda x: x.reshape(5,1), vs2)
  arr = np.concatenate(tuple(vs3), axis=1)
  per = np.percentile(loadings, 97)
  exceeds = filter(lambda x: np.max(x[1]) > per, enumerate(arr))
  exceeds_inds = map(lambda x: x[0], exceeds)
  exceeds_hgncs = map(lambda x: hgncs[x], exceeds_inds)
