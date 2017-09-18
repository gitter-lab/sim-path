import itertools as it
import os
import math
import numpy as np
import networkx as nx
import sys

def pathways_to_mat(pathways, id_to_index):
  """ Represent gene pathway membership as a matrix with shape (n_genes, n_pathways)

  Parameters
  ----------
  pathways : list of nx.Graph
    pathways used in simulation

  id_to_index : dict
    mapping of gene identifier to a numpy array index

  Returns
  -------
  mat : np.array
    matrix where each column is associated with a pathway in <pathways> and has a 1 in each
    row where that gene is in the pathway

  TODO
  ----
  For a single pathway, this representation has n_genes_in_pathway "units of signal" while
  the sampled signal has units <= n_genes_in_pathway (because a subset is sampled)
  """
  n_pathways = len(pathways)
  mat = np.zeros((len(id_to_index), n_pathways))
  for pathway_ind in range(len(pathways)):
    pathway = pathways[pathway_ind]
    for node in pathway.nodes_iter():
      index = id_to_index.get(node)
      if index is None:
        sys.stderr.write("Unrecognized gene ID: {}\n".format(node))
      else:
        mat[index,pathway_ind] = 1
  return mat

def argmax_nd(mat):
  """
  Return a tuple of the same length as the shape of mat which contains indices of argmax
  """
  argmax = np.argmax(mat)
  return np.unravel_index(argmax, mat.shape)

def remove_row_and_col(mat, row, col):
  mat = np.delete(mat, row, axis=0)
  mat = np.delete(mat, col, axis=1)
  return mat

def mat_to_bipartite(mat):
  G = nx.Graph()
  row_node_ids = []
  col_node_ids = []
  n_rows, n_cols = mat.shape
  for i in range(n_rows):
    node_i = "r{}".format(i)
    row_node_ids.append(node_i)
    G.add_node(node_i, bipartite=0)

  for j in range(n_cols):
    node_j = "c{}".format(j)
    col_node_ids.append(node_j)
    G.add_node(node_j, bipartite=1)

  for i in range(n_rows):
    for j in range(n_cols):
      G.add_edge("r{}".format(i), "c{}".format(j), {'weight': mat[i,j]})
  return G, row_node_ids, col_node_ids

def to_unit_vector(vec):
  rv = None
  norm = np.linalg.norm(vec)
  if(norm == 0):
    rv = vec
  else:
    rv = vec / norm
  return rv

def match(W_mat, pathways_mat):
  """
  Match latent factors to pathways

  Parameters
  ----------
  W_mat : np.array

  pathways_mat : np.array

  Returns
  -------
  rv : list of tpl
    each tpl is of the form (<latent_factor_id>, <pathway_id>, <distance>)
    where W_mat latent factors are identified with names "r0", "r1", etc. (for "row" of the 
    distance matrix constructed in this function) associated with each latent factor in W_mat; 
    and pathways are identified with names "c0", "c1", etc. (for "column" of the distance 
    matrix)

  TODO
  ----
  other type of transformation of distance to similarity? RBF?
  implement version that uses hypergeometic p-values as distance measure?
  """
  n_genes, n_factors = W_mat.shape
  n_genes2, n_pathways = pathways_mat.shape
  if(n_genes != n_genes2):
    raise ValueError("n_genes != n_genes2 : {} != {}".format(n_genes, n_genes2))
  if(n_factors != n_pathways):
    raise ValueError("n_factors != n_pathways: {} != {}".format(n_factors, n_pathways))

  dists = np.zeros((n_factors, n_pathways))
  max_unit_vec_dist = math.sqrt(2)
  for i in range(n_factors):
    factor_vec = W_mat[:,i]
    # scale distance by size of latent factor and pathway
    factor_vec_unit = to_unit_vector(factor_vec)
    for j in range(n_pathways):
      pathway_vec = pathways_mat[:,j]
      pathway_vec_unit = to_unit_vector(pathway_vec)
      dists[i,j] = np.linalg.norm(factor_vec_unit - pathway_vec_unit)
  dists_inv = 1 - (dists / max_unit_vec_dist)

  # use networkx for max weight matching
  G, factor_node_ids, pathway_node_ids = mat_to_bipartite(dists_inv)
  matching = nx.max_weight_matching(G)

  # I don't like the return value from networkx, reorganize data
  rv = []
  for factor_id in factor_node_ids:
    # not all nodes used in matching necessarily
    if factor_id in matching:
      pathway_id = matching[factor_id]
      dist = G[factor_id][pathway_id]['weight']
      tpl = (factor_id, pathway_id, dist)
      rv.append(tpl)
    else:
      # TODO temporary: why are they not included?
      weights = []
      for neighbor in G.neighbors(factor_id):
        weights.append(G[factor_id][neighbor]['weight'])
      sys.stderr.write("Factor {} not included in matching; weights: ".format(factor_id) + " ; ".join(map(str, weights)) + "\n")

  return rv

def parse_index_map(index_to_id_fp):
  id_to_index = {}
  index_to_id = []
  line_no = 0
  with open(index_to_id_fp) as fh:
    for line in fh:
      line = line.rstrip()
      index_to_id.append(line)
      id_to_index[line] = line_no
      line_no += 1
  return id_to_index

def normalize_num_pathways(W_mat, pathways_mat):
  """
  Prior to matching latent factors to pathways, this function MAY be invoked to transform
  W_mat or pathways_mat (whichever is smaller) to make it so the number of latent factors
  is the same as the number of pathways by adding a binary vector of all zeros (the empty pathway).

  Parameters
  ----------
  W_mat : np.array
    array with latent factors as columns and genes as rows with gene loadings into latent factors in cells

  pathways_mat : np.array
    array with pathways as columns and genes as rows with a 1 or 0 in each cell indicating whether that
    gene is present in the pathway or not

  Returns
  -------
  W_mat : np.array
    See description

  pathways_mat : np.array
    See description
  """
  # reshape 1D vectors if needed
  shp = W_mat.shape
  if(len(shp) == 1):
    W_mat = W_mat.reshape((shp[0], 1))
  shp = pathways_mat.shape
  if(len(shp) == 1):
    pathways_mat = pathways_mat.reshape((shp[0], 1))

  n_genes, n_latent_factors = W_mat.shape
  shp = pathways_mat.shape
  n_genes2, n_pathways = pathways_mat.shape
  if(n_genes != n_genes2):
    raise ValueError("{} != {}: number of genes identified by W_mat is not the same as those identified by pathways_mat".format(n_genes, n_genes2))

  arrs = [W_mat, pathways_mat]
  small_arr = -1
  large_arr = -1
  diff = 0
  if(n_latent_factors < n_pathways):
    small_arr = 0
    large_arr = 1
    diff = n_pathways - n_latent_factors
  elif(n_pathways < n_latent_factors):
    small_arr = 1
    large_arr = 0
    diff = n_latent_factors - n_pathways
  # else do nothing

  cols = np.zeros((n_genes, diff))
  arrs[small_arr] = np.concatenate((arrs[small_arr], cols), axis=1)

  return tuple(arrs)

def pathway_distance(pathways):
  """
  Compute distance between pairs of pathways

  Parameters
  ----------
  pathways : list of networkx.Graph
    each graph must have a graph attribute called "name" which is the name of the pathway 
    and will be used as one node class

  Returns
  -------
  dist_mat : numpy matrix with shape (n_pathways, n_pathways) 
    distance from one pathway to another according to an inverted Jaccard coefficient: 1 minus the number of
    genes shared in both pathways divided by the number of genes in either pathway

  pathway_nodes : list of str
    labels for either axis of dist_mat, in alphabetical sorted order
  """
  G = nx.Graph()
  pathway_nodes = []
  gene_nodes = []
  for pathway in pathways:
    pathway_name = pathway.graph['name']
    G.add_node(pathway_name, bipartite=0)
    for node in pathway:
      G.add_node(node, bipartite=1)
      G.add_edge(pathway_name, node)

  for node in G.nodes_iter():
    if G.node[node]['bipartite'] == 0:
      pathway_nodes.append(node)
    else:
      gene_nodes.append(node)
  pathway_nodes = sorted(pathway_nodes)

  pathway_node_ind_iter = it.combinations(range(len(pathway_nodes)), 2)
  dist_mat = np.zeros((len(pathway_nodes), len(pathway_nodes)))
  for i1, i2 in pathway_node_ind_iter:
    p1 = pathway_nodes[i1]
    p2 = pathway_nodes[i2]
    # jaccard coefficient function takes an iterator over pairs of node IDs but we
    # would like to refer to these nodes instead by their index in pathway_nodes for the
    # purpose of constructing dist_mat; this inner loop only has one iteration
    nx_iter = nx.jaccard_coefficient(G, [(p1, p2)])
    for p1, p2, jaccard_coef in nx_iter:
      dist_mat[i1,i2] = 1 - jaccard_coef
    
  return dist_mat, pathway_nodes

def get_pathways_in_dir(pathways_dir):
  """
  Returns
  -------
  pathways : list of networkx.Graph
  """
  fnames = os.listdir(pathways_dir)
  pathways = []
  for fname in fnames:
    basename, ext = os.path.splitext(fname)
    fpath = os.path.join(pathways_dir, fname)
    if ext == ".graphml":
      pathway = nx.read_graphml(fpath)
      pathway.graph['name'] = basename
      pathways.append(pathway)
  return pathways
