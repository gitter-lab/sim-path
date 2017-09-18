#!/usr/bin/env python
import argparse
import numpy as np
import networkx as nx
import sys
import sklearn.metrics
from ppi import simulation as sim
from ppi import ensembl
from ppi import parsers
from ppi import string_db
from ppi import kegg

def matching_id_to_ind(factor_id):
  return int(factor_id[1:])

def kegg_pathway_node_type(pathway):
  """
  Parameters
  ----------
  pathway : nx.Graph
    KEGG pathway downloaded by kegg.R
    KEGG's "hsa" identifiers have been converted (where possible, some hsa remain) to some
    other node type which is returned by this method

  Returns
  -------
  id_type : "hgnc_symbol", None
    if None, node identifiers are not recognized
  """
  # TODO
  id_type = None
  return id_type

# TODO move this to lib
def detect_and_parse_ppa_db(fh):
  db_type = parsers.detect_db(fh)
  G_ppa = None
  if(db_type == 'string'):
    G_ppa = string_db.parse_string_fh(fh)
  elif(db_type == 'irefindex'):
    G_ppa = iref.get_alt_id_ppi(fh)
  else:
    raise ArgumentError("Unrecognized protein database")
  return G_ppa, db_type

def main():
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("--index-to-id", "-i", help="Newline delimited file where each line contains one word, a gene identifier, which acts as a row label for <--W-path>.")
  parser.add_argument("--ppi-db", "-d", help="Either STRINGdb or iRefIndex database file", type=argparse.FileType('r'))
  parser.add_argument("--hgnc-to-ensp-map", "-m", help="Tab-separated value file with 4 fields: hgnc, ensg, enst, ensp", type=argparse.FileType('r'))
  parser.add_argument("--W-path", "-W", help="CSV file containing latent factor loadings for each gene", required=True)
  parser.add_argument("--pathways", "-p", nargs="+", help="KEGG graphml pathways used for simulation")
  parser.add_argument("--pathways-file", "-f", help="KEGG graphml pathways used for simulation provided in a newline-delimited file")
  args = parser.parse_args()

  pathway_fps = None
  if(args.pathways is None and args.pathways_file is None):
    sys.stderr.write("[error] One of --pathways or --pathways-file is required but neither is provided\n")
    sys.exit(21)
  elif(args.pathways is None and args.pathways_file is not None):
    pathway_fps = []
    with open(args.pathways_file, 'r') as fh:
      for line in fh:
        line = line.rstrip()
        pathway_fps.append(line)
  elif(args.pathways is not None and args.pathways_file is None):
    pathway_fps = args.pathways
  else:
    sys.stderr.write("[error] Exactly one of --pathways or --pathways-file is required but both are provided\n")
    sys.exit(22)
    
  id_to_index = sim.parse_index_map(args.index_to_id)
  W_mat = np.genfromtxt(args.W_path, delimiter=",")

  pathways = map(nx.read_graphml, pathway_fps)

  # transform pathway ids if needed
  if args.hgnc_to_ensp_map is not None and args.ppi_db is not None:
    # assume pathways have mostly HGNC symbols as identifiers
    G_ppa, db_type = detect_and_parse_ppa_db(args.ppi_db)
    if(db_type != 'string'):
      # this suggests a logical error, fail
      sys.stderr.write("[error] db_type is {} and not 'string'\n".format(db_type))
      sys.exit(23)
    hgnc_to_ensps_map = ensembl.map_hgnc_to_ensps(args.hgnc_to_ensp_map)
    for i in range(len(pathways)):
      pathway = pathways[i]
      mapped_obj = kegg.map_hgnc_to_ensps_graph(G_ppa, pathway, hgnc_to_ensps_map)
      pathways[i] = mapped_obj['G_sub']

  # need flag but also default would be nice
  pathways_mat = sim.pathways_to_mat(pathways, id_to_index)
  W_mat, pathways_mat = sim.normalize_num_pathways(W_mat, pathways_mat)
  matching = sim.match(W_mat, pathways_mat)

  # 1) evaluate performance by a distance metric between latent factor and pathway
  # dists are in range 0 to 1
  print("\t".join(["#" + "factor_id", "pathway_id", "distance"]))
  total_dist = 0
  n_dists = 0
  for match in matching:
    factor_id, pathway_id, dist = match
    print("\t".join(map(str, match)))
    total_dist += dist
    n_dists += 1
  print(total_dist / n_dists)

  # 2) evaluate performance by ROC
  print("\t".join(["#" + "factor_id", "pathway_id", "auc"]))
  total_auc = 0
  n_auc = 0
  for match in matching:
    factor_id, pathway_id, dist = match
    factor_id = matching_id_to_ind(factor_id)
    pathway_id = matching_id_to_ind(pathway_id)

    y_score = W_mat[:,factor_id]
    y_true = pathways_mat[:,pathway_id]
    precision, recall, thresholds = sklearn.metrics.precision_recall_curve(y_true, y_score)
    auc = sklearn.metrics.average_precision_score(y_true, y_score)
    print("\t".join(map(str, [factor_id, pathway_id, auc])))
    total_auc += auc
    n_auc += 1
  print(total_auc / n_auc)

if __name__ == "__main__":
  main()
