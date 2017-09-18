#!/usr/bin/env python
import argparse, sys
import networkx as nx
from ppi import irefindex as iref
from ppi import string_db as sdb
from ppi import parsers
from ppi import hgnc
from ppi import kegg
from ppi import ensembl
import numpy.random as nprand
import os, os.path

def subtract_pathways(G_ppi, pathways):
  """
  Remove <pathways> from <G_ppi>

  Parameters
  ----------
  G_ppi : nx.Graph
    NOTE this graph will be mutated

  pathways: list of nx.Graph

  Returns
  -------
  G_minus : nx.Graph
  """
  edges = set() # set of 2-tuples
  for pathway in pathways:
    for edge in pathway.edges_iter():
      edges.add(edge)
    
  G_ppi.remove_edges_from(edges)
  G_minus = G_ppi # TODO dont mutate G_ppi?
  return G_minus

def sample_gene_list(pathways, list_size, seed=None):
  """
  Sample from a multinomial distribution the number of genes from each pathway that will be used
  Then, from each pathway, sample that number of genes without replacement.

  Parameters
  ----------
  pathways : list of nx.Graph
    last graph may be a noise graph: e.g. setting to global PPI minus all pathways

  Returns
  -------
  gene_list : list of str
  """
  gene_list = []
  sample_sizes = nprand.multinomial(list_size, [1/float(len(pathways))] * len(pathways))
  for i in range(len(pathways)):
    pathway = pathways[i]
    # TODO this may choose nodes which failed conversion
    # TODO alternative to trying to sample more nodes than are in the pathway?
    sample = []
    try:
      sample = nprand.choice(pathway.nodes(), size=sample_sizes[i], replace=False)
    except ValueError as err:
      # then we tried to sample more nodes than are in the pathway
      sample = pathway.nodes()
    gene_list += list(sample)
  return list(set(gene_list))

def main():
  parser = argparse.ArgumentParser(description="""
Simluate genetic screens from KEGG pathways and iRefIndex
""")
  # TODO noise currently not used
  parser.add_argument("--noise", "-n", type=float, help="Value in [0,1) providing the probability that a node is sampled from an iRefIndex vertex that is not in one of the pathways", default=0)
  parser.add_argument("--list-size", "-s", type=int, help="Number of nodes sampled for each pathway", default=300)
  parser.add_argument("--m-gene-lists", "-m", type=int, help="Number of gene lists to create", default=10)
  parser.add_argument("--outdir", "-o", help="Directory to write simulated genetic screens to", default=os.curdir)
  parser.add_argument("--kegg-pathways", help="graphml files of KEGG pathways with identifiers as HGNC gene symbols", nargs='+')
  parser.add_argument("--pathways-file", help="file containing a filepath of a KEGG pathway graphml file on each line")
  parser.add_argument("--ppi-db", help="iRefIndex or STRING database file", type=argparse.FileType('r'), required=True)
  parser.add_argument("--hgnc-ensp-map", help="Mapping file among HGNC, ENSG, ENST, and ENSP", type=argparse.FileType('r'))
  args = parser.parse_args()

  G_ppi = None
  hgnc_ensp_map = None
  db_type = parsers.detect_db(args.ppi_db)
  if(db_type == 'string' and args.hgnc_ensp_map is None):
    sys.stderr.write("Database type string requires mapping file --hgnc-ensp-map\n")
    sys.exit(3)
  if(db_type == 'irefindex'):
    G_ppi = iref.get_alt_id_ppi(args.ppi_db)
  elif(db_type == 'string'):
    G_ppi = sdb.parse_string_fh(args.ppi_db)
    hgnc_ensp_map = ensembl.map_hgnc_to_ensps(args.hgnc_ensp_map)
  # otherwise error raised

  pathways_files = []
  if(args.kegg_pathways is None and args.pathways_file is None):
    sys.stderr.write("[error] Either --kegg-pathways or --pathways-file must be provided and neither were provided\n")
    sys.exit(1)
  elif(args.kegg_pathways is None and not args.pathways_file is None):
    with open(args.pathways_file) as fh:
      for line in fh:
        line = line.rstrip()
        pathways_files.append(line)
  elif(not args.kegg_pathways is None and args.pathways_file is None):
    pathways_files = args.kegg_pathways
  else:
    sys.stderr.write("[error] Cannot provide both --kegg-pathways and --pathways-file\n")
    sys.exit(2)
  pathways = map(lambda x: nx.read_graphml(x), pathways_files)
  if(db_type == 'string'):
    # Lossy translate pathway identifiers to Ensembl
    for i in range(len(pathways)):
      pathway = pathways[i]
      map_obj = kegg.map_hgnc_to_ensps_graph(G_ppi, pathway, hgnc_ensp_map)
      pathways[i] = map_obj['G_sub']

  G_minus = subtract_pathways(G_ppi, pathways)
  pathways_and_noise = list(pathways) + [G_minus]

  gene_lists = []
  for i in range(args.m_gene_lists):
    gene_list = sample_gene_list(pathways_and_noise, args.list_size)
    gene_lists.append(gene_list)

  list_no = 1
  for gene_list in gene_lists:
    list_name = "list{}.txt".format(list_no)
    fp = os.path.join(args.outdir, list_name)
    with open(fp, 'w') as fh:
      for gene in gene_list:
        fh.write("{}\n".format(gene))
    list_no += 1

if __name__ == "__main__":
  main()
