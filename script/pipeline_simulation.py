#!/usr/bin/env python
import numpy as np
import sys, argparse
import os, os.path
import subprocess
from ppi import script_utils
import networkx as nx

def uniform_random_pathways(pathways_dir, n_pathways, seed=None):
  """
  Returns
  -------
  sampled_pathways : list of str
    list of filepaths to pathways to use for the simulation
  """
  if seed is not None:
    np.random.seed(seed)
  pathways = os.listdir(pathways_dir)
  pathways = list(map(lambda x: os.path.join(pathways_dir, x), pathways))
  sampled_pathways = list(np.random.choice(pathways, n_pathways))
  return sampled_pathways

def distance_random_pathways(pathways_dir, n_pathways, pathways_file, job_graph, job_id, outdir, seed=None):
  """
  Returns
  -------
  sampled_pathways : list of str
    list of filepaths to pathways to use for the simulation

  See also uniform_random_pathways
  """
  # select uniformly at random one pathway
  # then select n_pathways - 1 more pathways based on their distance to the selected pathway
  attrs = {
    'exe': "simulate_pathways.py",
    'args': ["--n-pathways", str(n_pathways), "-o", pathways_file, pathways_dir],
    'out': os.path.join(outdir, "simulate_pathways.out"),
    'err': os.path.join(outdir, "simulate_pathways.err")
  }
  if seed is not None:
    attrs['args'] += ["--seed", str(seed)]
  job_graph.add_node(job_id, attrs)

def main():
  sys.stdout.write(" ".join(sys.argv) + "\n")
  parser = argparse.ArgumentParser(description="Simulation pipeline")
  # arguments are grouped according to their influence on a stage of the pipeline
  parser.add_argument("--n-simulations", type=int, help="Number of simulations to run", default=1)

  parser.add_argument("--n-pathways", type=int, help="Number of pathways to sample in each simulation", default=3)
  parser.add_argument("--pathways-range", "-p", type=int, help="For example, if the number of <kegg_pathways> provided is 3 and <--pathways-range> is 2, then factorizations with k in {1,2,3,4,5} will be performed", default=0)
  parser.add_argument("--pathway-method", type=str, default="distance", help="Procedure for selecting which pathways to simulate genetic screen hits from. Must be either \"distance\" or \"uniform\".")
  parser.add_argument("--threshold", type=float, default=0.0, help="Confidence threshold for edges in the PPA network. Must be in [0,1]")
  # TODO pass through options for simulate_pathways
  #parser.add_argument("--distance-min", "-d", type=float, default=0.0)
  #parser.add_argument("--distance-max", "-D", type=float, default=1.0)
  parser.add_argument("--hgnc-ensp-map", type=str, help="Mapping file among HGNC, ENSG, ENST, and ENSP")

  parser.add_argument("--n-gene-lists", "-g", type=int, default=10)
  parser.add_argument("--gene-list-range", "-r", type=int, default=0, help="See pathways-range; same idea but for the number of gene lists")
  parser.add_argument("--same-lists", "-s", action="store_true", help="If this flag is provided, varying numbers of the same gene lists will be used rather than simulating new gene lists for each run")

  # misc options
  parser.add_argument("--seed", type=int, default=None, help="Seed for random number generator")
  parser.add_argument("--dry-run", action="store_true", help="If provided, prepare DAG but do not run it")
  # TODO check mutex local, chtc, missing
  parser.add_argument("--local", action="store_true", help="If provided, run pipeline locally instead of on Condor")
  parser.add_argument("--chtc", action="store_true", help="If provided, prepare files for Condor by including transfer_input_files in the Condor submit file")

  parser.add_argument("ppi_db", type=str, help="iRefIndex or STRING database file")
  parser.add_argument("kegg_pathways_dir")
  parser.add_argument("outdir", type=str)
  args = parser.parse_args()

  # specify command prerequisites
  # store arguments to script_utils.run_command as a node attribute
  job_graph = nx.DiGraph()
  job_id = 0
  for sim_i in range(args.n_simulations):
    sim_outdir = os.path.join(args.outdir, "s{}".format(sim_i))
    script_utils.mkdir_p(sim_outdir)

    # sample pathways from all KEGG pathways: prepare pathways_file for simulate_screens.py
    pathways_file = os.path.join(sim_outdir, "simulation_pathways.txt")
    pathways_job_id = None
    if(args.pathway_method == "distance"):
      # TODO additional arguments to distance_random_pathways
      distance_random_pathways(args.kegg_pathways_dir, args.n_pathways, pathways_file, job_graph, job_id, sim_outdir, seed=args.seed)
      pathways_job_id = job_id
      job_id += 1

    elif(args.pathway_method == "uniform"):
      kegg_pathways = uniform_random_pathways(args.kegg_pathways_dir, args.n_pathways, seed=args.seed)
      # log pathways used
      with open(pathways_file, 'w') as fh:
        for pathway in kegg_pathways:
          fh.write(pathway + "\n")
    else:
      sys.stderr.write("[error] pathway_method \"{}\" must be \"distance\" or \"uniform\"\n")
      sys.exit(1)

    # ground truth number of pathways
    k_pathways = args.n_pathways

    # simulate screens
    screen_fps = []
    simulate_job_id = None
    if(args.same_lists):
      # then simulate screens once for use by the rest of the pipeline
      attrs = {
        'exe': "simulate_screens.py",
        'args': ["-o", sim_outdir, "-m", str(args.n_gene_lists + args.gene_list_range), "--ppi-db", args.ppi_db, "--pathways-file", pathways_file],
        'out': os.path.join(sim_outdir, "simulate_screens.out"),
        'err': os.path.join(sim_outdir, "simulate_screens.err")
      }
      if(args.hgnc_ensp_map is not None):
        attrs['args'] += ["--hgnc-ensp-map", args.hgnc_ensp_map]
      job_graph.add_node(job_id, attrs)
      simulate_job_id = job_id
      job_id += 1
      if(pathways_job_id is not None):
        job_graph.add_edge(pathways_job_id, simulate_job_id)
      for i in range(1, args.n_gene_lists + args.gene_list_range + 1):
        screen_fps.append(os.path.join(sim_outdir, "list{}.txt".format(i)))

    # run nbs
    # use various numbers of gene lists
    min_g = args.n_gene_lists - args.gene_list_range if (args.n_gene_lists - args.gene_list_range > 0) else 1
    max_g = args.n_gene_lists + args.gene_list_range
    for g in range(min_g, max_g+1, 1):
      # determine gene lists
      g_outdir = os.path.join(sim_outdir, "g{}".format(g))
      script_utils.mkdir_p(g_outdir)
      this_screen_fps = []
      if(args.same_lists):
        # then use a subset of the already generated gene lists
        this_screen_fps = screen_fps[0:g]
      else:
        # then simulate screens for each run of the pipeline
        this_screen_fps = []
        attrs = {
          'exe': "simulate_screens.py", 
          'args': ["-o", g_outdir, "-m", str(g), args.ppi_db] + kegg_pathways,
          'out': os.path.join(g_outdir, "simulate_screens.out"),
          'err': os.path.join(g_outdir, "simulate_screens.err")
        }
        job_graph.add_node(job_id, attrs)
        simulate_job_id = job_id
        job_id += 1
        if(pathways_job_id is not None):
          job_graph.add_edge(pathways_job_id, simulate_job_id)
        for i in range(1, g+1, 1):
          this_screen_fps.append(os.path.join(g_outdir, "list{}.txt".format(i)))

      # use gene lists for diffusion
      diffused_path = os.path.join(g_outdir, "diffused.mat")
      ids_path = os.path.join(g_outdir, "diffused_ids.txt")
      attrs = {
        'exe': "nbs.py",
        'args': [
            "--ids", ids_path, 
            "--mat-file", diffused_path,
            "--threshold", str(args.threshold)
          ] + [args.ppi_db] + this_screen_fps,
        'out': os.path.join(g_outdir, "nbs.out"),
        'err': os.path.join(g_outdir, "nbs.err")
      }
      job_graph.add_node(job_id, attrs)
      nbs_job_id = job_id
      job_graph.add_edge(simulate_job_id, nbs_job_id)
      job_id += 1

      # run gnmf
      # test various numbers of k_pathways, including the ground truth number
      # TODO this leaves behind matlab runtime directory in location where this script was run from
      # TODO wu_init flag
      min_k = k_pathways - args.pathways_range if (k_pathways - args.pathways_range > 0) else 1
      max_k = k_pathways + args.pathways_range
      for k in range(min_k, max_k+1, 1):
        k_outdir = os.path.join(g_outdir, "k{}".format(k))
        script_utils.mkdir_p(k_outdir)
        W_path = os.path.join(k_outdir, "W.csv")
        H_path = os.path.join(k_outdir, "H.csv")
        res_path = os.path.join(k_outdir, "residual.csv")
        runtime_path = os.path.join(os.environ['SIMPATH_REPO_DIR'], "bin", "r2014b.tar.gz")
        gnmf_script_path = os.path.join(os.environ['SIMPATH_REPO_DIR'], "bin", "gnmf_script")
        attrs = {
          'exe': "mcc_script.sh",
          'args': ["-m", runtime_path, "--", gnmf_script_path, diffused_path, W_path, H_path, res_path, str(k)],
          'out': os.path.join(k_outdir, "run_gnmf_script.out"),
          'err': os.path.join(k_outdir, "run_gnmf_script.err"),
          'requirements': 'OpSysMajorVer == 7'
        }
        gnmf_job_id = job_id
        job_graph.add_node(gnmf_job_id, attrs)
        job_graph.add_edge(nbs_job_id, gnmf_job_id)
        job_id += 1

        # run evaluation
        eval_job_id = job_id
        
        attrs = {
          'exe': "evaluate_simulation.py",
          'args': ["-i", ids_path, "-d", args.ppi_db, "-m", args.hgnc_ensp_map, "-W", W_path, "-f", pathways_file],
          'out': os.path.join(k_outdir, "evaluate_simulation.out"),
          'err': os.path.join(k_outdir, "evaluate_simulation.err")
        }
        job_graph.add_node(eval_job_id, attrs)
        job_graph.add_edge(gnmf_job_id, eval_job_id)
        job_id += 1

  # end simulation loop
  # execute job graph
  condor = (not args.local)
  job_ids = script_utils.run_digraph(args.outdir, job_graph, condor=condor, dry_run=args.dry_run)
  for job_id in job_ids:
    print(job_id)


if __name__ == '__main__':
  main()
