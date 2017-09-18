#!/usr/bin/env python
import numpy as np
import sys, argparse
import random
import os, os.path
from ppi import simulation as sim

# TODO what if two "far" pathways are "close" to each other
# need something like max weight clique which is NP-complete
# but complete graph on KEGG pathways as nodes may be sufficiently small

def main():
  parser = argparse.ArgumentParser(description="""
Choose a set of pathways in <pathways_dir> based on their proximity.
""")
  parser.add_argument("--seed", default=None, type=int, help="Seed for random selection if <--source-pathway> is not provided")
  parser.add_argument("--source-pathway", "-s", default=None, type=str, help="Pathway to calculate distances from. If not provided, one will be selected uniformly at random from the pathways_dir")
  parser.add_argument("--n-pathways", "-n", default=5, type=int, help="Number of pathways to sample")
  parser.add_argument("--min-distance", "-d", default=0, type=float, help="Minimum distance from <--source-pathway> to be eligible for selection. Expressed as a percentile of all pathway distances from the source pathway. Inclusive. Default=0")
  parser.add_argument("--max-distance", "-D", default=1, type=float, help="Maximum distance from <--source-pathway> to be eligible for selection. Expressed as a percentile of all pathway distances from the source pathway. Inclusive. Default=1")
  parser.add_argument("--increasing-order", "-i", action='store_true', help="If provided, reverse the default decreasing order selection in range.")
  parser.add_argument("--outfile", "-o", default=sys.stdout, type=argparse.FileType("w"))
  parser.add_argument("pathways_dir")
  args = parser.parse_args()

  pathways = sim.get_pathways_in_dir(args.pathways_dir)
  # TODO this currently computes the distance between all pairs of pathways which is unnecessary here
  # but is something I plan to use in the future for example the clique TODO above
  dist_mat, pathway_nodes = sim.pathway_distance(pathways)
  source_pathway_index = None
  if(args.source_pathway is None):
    source_pathway_index = random.randint(0, len(pathway_nodes)-1)
  else:
    # then args.source_pathway is the pathway name for a file in pathways_dir
    try:
      source_pathway_index = pathway_nodes.index(args.source_pathway)
    except ValueError as err:
      sys.stderr.write("Invalid source pathway \"{}\". Must be a name of a file in pathways_dir \"{}\" (without extension).\n".format(args.source_pathway, args.pathways_dir))
      sys.exit(1)

  dists = dist_mat[source_pathway_index]
  source_pathway_id = pathway_nodes[source_pathway_index]

  reverse = True
  if(args.increasing_order):
    reverse = False
  dist_and_id = sorted(zip(dists, pathway_nodes, range(len(dists))), key=lambda x: x[0], reverse=reverse)

  def filt_fn(x):
    min_v = np.round(np.percentile(dists, args.min_distance*100), decimals=3)
    max_v = np.round(np.percentile(dists, args.max_distance*100), decimals=3)
    x_v = np.round(x[0], decimals=3)
    if(min_v <= x_v and x_v <= max_v):
      return True
    else:
      return False
  dist_and_id_filt = filter(filt_fn, dist_and_id)

  # TODO this block assumes that one of the pathways in _filt is not the source pathway
  select_dist_and_id = None
  if(len(dist_and_id_filt) < args.n_pathways-1):
    # then warn, but guarantee n_pathways TODO
    sys.stderr.write("[warning] insufficient pathways ({} of {}) in range [{}, {}]\n".format(len(dist_and_id_filt), args.n_pathways, args.min_distance, args.max_distance))
    diff = args.n_pathways - len(dist_and_id_filt)
    select_dist_and_id = dist_and_id_filt[:args.n_pathways-1]
  else:
    select_dist_and_id = dist_and_id_filt[:args.n_pathways-1]

  pathway_ids = list(map(lambda x: x[1], select_dist_and_id))
  if(source_pathway_id not in pathway_ids):
    # if looking for closest pathways, then tie for 0 dist
    # otherwise, this case is expected, remove furthest to make room for source
    if(reverse):
      pathway_ids = pathway_ids + [source_pathway_id]
    else:
      pathway_ids = [source_pathway_id] + pathway_ids

  # TODO what if not graphml
  for pathway_id in pathway_ids:
    args.outfile.write(os.path.join(args.pathways_dir, pathway_id + ".graphml") + "\n")

if __name__ == "__main__":
  main()
