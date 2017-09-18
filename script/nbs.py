#!/usr/bin/env python
import scipy.io
import numpy as np
import argparse, sys
from ppi import irefindex as iref
from ppi import parsers
from ppi import string_db as sdb
import stratipy.filtering_diffusion as diffuse
import scipy.sparse as sp
import sys
import os.path
import pickle

def parse_gene_list_file(fh):
  rv = []
  for line in fh:
    line = line.rstrip()
    rv.append(line)
  return rv

def check_or_set_fh(arg, default_base):
  """Convenience function for setting default file handle
  """
  rv = None
  if(arg is None):
    path = os.path.join(os.curdir, default_base)
    if(os.path.exists(path)):
      sys.stderr.write("Refusing to overwrite {}\n".format(path))
      exit(2)
    else:
      rv = open(path, "w")
  else:
    rv = arg
  return rv

def main():
  # TODO currently only does smoothing step
  parser = argparse.ArgumentParser(description="Run network based stratification")
  parser.add_argument("ppi_db", type=str, help="irefindex database file")
  parser.add_argument("gene_lists", nargs="+", help="one or more files with an HGNC identifier on each line")
  parser.add_argument("--alpha", "-a", type=float, default=0.7,
    help="Diffusion rate parameter")
  parser.add_argument("--tolerance", "-t", type=float, default=10e-6,
    help="Tolerance threshold for diffusion; stop when change in diffused matrix crosses below threshold")
  parser.add_argument("--threshold", type=float, default=0.0),
  parser.add_argument("--ids", type=argparse.FileType("w"),
    help="File to write association between an index in the smoothed data matrix and a gene identifier")
  parser.add_argument("--mat-file", "-m", nargs=1, type=argparse.FileType("w"),
    help=".mat file to write 'smoothed_mat' and 'laplacian' matrices to; if provided, laplacian and diffused are ignored")
  parser.add_argument("--laplacian", "-K", type=argparse.FileType("w"),
    help="File to write graph Laplacian to")
  parser.add_argument("--diffused", "-d", type=argparse.FileType("w"),
    help="Diffused matrix")
  args = parser.parse_args()

  # set default output files but dont overwrite if a default is used (overwrite if file is explicitly mentioned)
  args.ids = check_or_set_fh(args.ids, "diffused_ids.txt")
  do_mat = False
  if(args.mat_file is None):
    args.laplacian = check_or_set_fh(args.laplacian, "laplacian.csv")
    args.diffused = check_or_set_fh(args.diffused, "smoothed_mat.csv")
  else:
    do_mat = True

  # (1) compute smoothed matrix
  gene_lists = []
  for gene_path in args.gene_lists:
    with open(gene_path) as fh:
      gene_lists.append(parse_gene_list_file(fh))

  fh = open(args.ppi_db)
  db_type = parsers.detect_db(fh)
  G_ppi = None
  if(db_type == 'string'):
    # TODO apply threshold only for laplacian?
    G_ppi = sdb.parse_string_fh(fh, threshold=args.threshold)
  elif(db_type == 'irefindex'):
    # TODO threshold
    G_ppi = iref.get_alt_id_ppi(fh) # defaults to node IDs as HGNC symbols
  # else error
  adj, ids = iref.get_adj_mat(G_ppi)
  fh.close()

  # bind first argument of get_row_vec in closure
  def get_row_vec_for_gene_list(gene_list):
    row_vec, missing = iref.get_row_vec_alt(ids, gene_list)
    sys.stderr.write("missing {}/{} HGNC identifiers: {}\n".format(len(missing), len(gene_list), ", ".join(missing)))
    return row_vec
  row_vecs = map(get_row_vec_for_gene_list, gene_lists)

  mat = sp.vstack(row_vecs)
  sys.stderr.write("mat.shape: {}\n".format(mat.shape))
  sys.stderr.write("adj.shape: {}\n".format(adj.shape))
  smoothed_mat = diffuse.propagation(mat, adj, alpha=args.alpha, tol=args.tolerance)

  # (2) compute laplacian 
  # TODO use a "influence distance" "kNN" graph instead of full graph laplacian
  laplacian, ids_lap = iref.get_laplacian(G_ppi)
  assert ids == ids_lap, "ids != ids_lap: smoothed matrix and Laplacian matrix do not agree on gene identities"

  # (3) write output files
  for id in ids:
    args.ids.write("{}\n".format(id))

  # take care in writing the laplacian; it is a n_gene x n_gene matrix (so approx 20k x 20k)
  # which is reasonably large; writing it as a .mat file is one option; another choice
  # is the MatrixMarket format which can write sparse matrices
  if(do_mat):
    # NOTE must verify that type is float because MatLab will interpret the result 
    # differently (incorrectly) if not
    scipy.io.savemat(args.mat_file[0], { 
      "smoothed_mat": smoothed_mat.astype(float),
      "laplacian": laplacian.astype(float)
    })
  else:
    scipy.io.mmwrite(args.laplacian, laplacian)
    scipy.io.mmwrite(args.diffused, smoothed_mat)

if __name__ == '__main__':
  main()
