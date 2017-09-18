"""
Functions for manipulating graphml output from KEGGgraph R library
"""
import networkx as nx
import itertools as it

def relabel_nodes(G_kegg):
  """
  The graphml file written by kegg.R uses consecutive non-negative integers as node identifiers
  and stores the gene identifier in node data. This function relabels nodes from their
  integer identifiers to their gene identifiers so that the identifiers match with iRefIndex.
  """
  mapping = {}
  for node in G_kegg.nodes_iter():
    node_data = G_kegg.node[node]
    mapping[node] = str(node_data['name'])

  G_kegg = nx.relabel_nodes(G_kegg, mapping)
  return G_kegg

def map_hgnc_to_ensps_graph(G_ppa, G_kegg, hgnc_to_ensps_map):
  """
  STRINGdb reports protein-protein assocations (interaction, phosphorylation, ubiquitination, etc.) between pairs of proteins.
  KEGG reports these assocations between pairs of genes.
  Gene -> Protein is a one to many mapping.
  For a gene association pair, few of the mapped protein pairs exist in the PPA database.
  Therefore, using the edges in G_ppa helps map G_kegg to more specific proteins.

  Parameters
  ----------
  G_ppa : nx.Graph

  G_kegg : nx.DiGraph

  hgnc_to_ensps_map : dict
    mapping of HGNC symbol to a set of ENSPs

  Returns
  -------
  G_sub : nx.Graph
  """
  G_sub = nx.Graph()
  ensps = set()
  missing_nodes = set()
  missing_edges = []
  for n0_hgnc, n1_hgnc in G_kegg.edges_iter():
    n0_ensps = hgnc_to_ensps_map.get(n0_hgnc)
    n1_ensps = hgnc_to_ensps_map.get(n1_hgnc)
    if(n0_ensps is None):
      missing_nodes.add(n0_hgnc)
    if(n1_ensps is None):
      missing_nodes.add(n1_hgnc)
    if(n0_ensps is not None and n1_ensps is not None):
      any_edge = False
      for n0_ensp, n1_ensp in it.product(n0_ensps, n1_ensps):
        if G_ppa.has_edge(n0_ensp, n1_ensp):
          ensps.add(n0_ensp)
          ensps.add(n1_ensp)
          weight = G_ppa[n0_ensp][n1_ensp]['weight']
          G_sub.add_edge(n0_ensp, n1_ensp, {'weight': weight})
          any_edge = True
      if not any_edge:
        missing_edges.append((n0_hgnc, n1_hgnc))
  return {
    'G_sub': G_sub,
    'ensps': ensps,
    'missing_nodes': missing_nodes,
    'missing_edges': missing_edges
  }
