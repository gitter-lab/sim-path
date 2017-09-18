"""
Functions for working with iRefIndex in terms of NetworkX graphs (instead of graph_tool)
"""

import numpy as np # for np.Inf TODO alternative?
import networkx as nx
import itertools
from .irefindex_parser import *

def get_ppi(fh, n_rows=None, hgnc=False):
  """Parse ppi network from irefindex

  Parameters
  ----------
  fh : file handle 
    irefindex data

  n_rows : int
    parse only first n_rows data rows from fh

  hgnc : boolean
    if true, add hgnc identifiers of interactors as node attributes

  Returns
  -------
  G : networkx undirected graph
    ppi network
  """

  try:
    fh.seek(0)
  except AttributeError:
    pass

  header = parse_header(fh)

  if n_rows is None:
    # then parse all data rows
    n_rows = np.Inf

  # Construct interaction graph
  G = nx.Graph()
  i = 0
  for line in fh:
    line = line.strip()
    datum = line.split("\t")
    if (not is_human_human(header, datum)):
      # include only human-human interactions
      continue

    node_a = None
    node_b = None
    id_a, id_b = get_checksums(header, datum)
    if (hgnc):
      # then add hgnc as a node attribute
      hgncs_a, hgncs_b = get_hgnc_pair(header, datum)
      attrs_a = { HGNCS_KEY: hgncs_a }
      attrs_b = { HGNCS_KEY: hgncs_b }
      G.add_node(id_a, attrs_a)
      G.add_node(id_b, attrs_b)

    G.add_edge(id_a, id_b)

    i += 1
    if i > n_rows:
      break
  return G

def get_alt_id_ppi(fh, id_type='hgnc', n_rows=None):
  """
  Construct a Protein-Protein interaction network from <fh> where nodes are identifiers of the type <id_type>.
  If an iRefIndex interactor identifier maps to multiple identifiers of the specified type, the resuling network
  contains all edges in the cartesian product of idsA and idsB.

  Parameters
  ----------
  fh : io-like
    irefindex database file

  id_type : string
    One of crogid, entrezgene/locuslink, genbank_protein_gi, hgnc, icrogid, irogid, pdb, refseq, rogid, uniprotkb
    see get_alias_names or KNOWN_ALIASES

  n_rows : int
    stop reading fh after <n_rows> lines

  Returns
  -------
  G : nx.Graph
  
  See also get_ppi_hgnc
  """
  try:
    fh.seek(0)
  except AttributeError:
    raise ArgumentError("fh must implement seek")

  if id_type not in KNOWN_ALIASES:
    raise ValueError("id_type: {} is invalid and not one of {}", id_type, ", ".join(KNOWN_ALIASES))

  min_lpr, max_lpr = get_lpr_extrema(fh)
  fh.seek(0)
  header = parse_header(fh)

  if n_rows is None:
    # then parse all data rows
    n_rows = np.Inf

  # Construct interaction graph
  G = nx.Graph()
  i = 0
  for line in fh:
    line = line.rstrip()
    datum = line.split("\t")
    lpr = get_lpr(header, datum)
    if lpr is None:
      # then assign worst confidence score
      lpr = max_lpr
    weight = lpr_to_conf(min_lpr, max_lpr, lpr)
    if (not is_human_human(header, datum)):
      # include only human-human interactions
      continue

    ids_a, ids_b = get_id_pair(header, datum, id_type=id_type)
    for id_pair in itertools.product(ids_a, ids_b):
      if G.has_edge(id_pair[0], id_pair[1]):
        # take the maximum edge weight ("confidence")
        edge_data = G.get_edge_data(id_pair[0], id_pair[1])
        other_weight = edge_data['weight']
        if weight > other_weight:
          # then update the edge weight, otherwise leave the edge as is
          G.add_edge(id_pair[0], id_pair[1], {'weight': weight})
      else:
        G.add_edge(id_pair[0], id_pair[1], {'weight': weight})

    i += 1
    if i > n_rows:
      break

  return G

def get_ppi_hgnc(fh, n_rows=None):
  """See get_ppi. As an alternative to hgnc=True, construct a network with
  the networkx node identifiers as hgnc symbols. This results in a different
  network than get_ppi(fh, hgnc=True) because not all the interactors in irefindex
  have hgnc identifiers.
  """
  try:
    fh.seek(0)
  except AttributeError:
    pass

  header = parse_header(fh)

  if n_rows is None:
    # then parse all data rows
    n_rows = np.Inf

  # Construct interaction graph
  G = nx.Graph()
  i = 0
  for line in fh:
    line = line.strip()
    datum = line.split("\t")
    if (not is_human_human(header, datum)):
      # include only human-human interactions
      continue

    hgncs_a, hgncs_b = get_hgnc_pair(header, datum)
    for hgnc_pair in itertools.product(hgncs_a, hgncs_b):
      G.add_edge(hgnc_pair[0], hgnc_pair[1])

    i += 1
    if i > n_rows:
      break

  return G

def components(fh, **ppi_opts):
  """Return a list of sets of connected components

  Parameters
  -----------
  fh : file handle
    irefindex data file

  ppi_opts : dict
    keyword arguments for get_ppi

  Returns
  -------
  comps_l : list of sets
    Each set member is a tuple for an iRefIndex interactor identifier (<id-type>, <id-value>) 
    e.g. ('rogid', 'NqjvDnObnt6a2DcQGJjMD2/mI3I9606')
  """
  G = get_ppi(fh, **ppi_opts)
  comps_g = nx.connected_components(G) # set generator
  comps_l = sorted(comps_g, key=len, reverse=True)
  return comps_l

def write_edge_per_line(G, ofh):
  """Write graph in form "<id_a> <id_b>\n"
  """
  for e in G.edges_iter():
    # nodes are tuples of (<id_type>, <id_value>)
    n1 = e[0]
    n2 = e[1]
    ofh.write("{}\t{}\n".format(n1[1], n2[1]))

def get_adj_mat(G):
  """Represent ppi network as adjacency matrix

  Parameters
  ----------
  G : networkx graph
    ppi network, see get_ppi()

  Returns
  -------
  adj : square sparse scipy matrix
    (i,j) has a 1 if there is an interaction reported by irefindex

  ids : list
    same length as adj, ith index contains irefindex unique identifier for gene whose interactions are
    reported in the ith row of adj
  """
  ids = G.nodes()
  adj = nx.to_scipy_sparse_matrix(G, nodelist=ids, dtype=bool)

  return adj, ids

def get_laplacian(G):
  """Get graph laplacian

  Parameters
  ----------
  G : networkx graph
    ppi network

  Returns
  -------
  laplacian : square scipy sparse matrix
    graph laplacian

  ids : list
    same length as lapacian containing node ids for each index
  """
  ids = G.nodes()
  laplacian = nx.laplacian_matrix(G, nodelist=ids)
  return laplacian, ids

def filter_node_ids(G, ids):
  """
  See irefindex_parser.filter_ids
  """
  found = []
  missing = []
  for id in ids:
    if id in G:
      found.append(id)
    else:
      missing.append(id)
  return (found, missing)
