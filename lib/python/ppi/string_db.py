"""
Functions to parse STRING protein-protein interaction database:
head -n2 9606.protein.links.full.v10.5.txt
protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
9606.ENSP00000000233 9606.ENSP00000263431 0 0 0 0 0 0 53 0 176 0 0 0 128 260
"""
import networkx as nx

def parse_line_abc(line):
  """
  Filter STRING data to just protein,protein,confidence and remove taxonomy code from protein identifiers
  """
  line = line.rstrip()
  words = line.split()
  if(len(words) != 16):
    raise RuntimeError("line {} in database file handle is not 16 tokens long".format(line_no))
  p1 = words[0]
  p2 = words[1]
  # remove human taxonomy code prefix
  p1 = trim_taxonomy_code(p1)
  p2 = trim_taxonomy_code(p2)
  score = words[-1]
  conf = int(score) / 1000.0
  return (p1, p2, conf)

def trim_taxonomy_code(protein_id):
  return protein_id[5:]

def parse_string_fh(fh, threshold=0.0):
  """
  Parameters
  ----------
  fh : file-like
    STRING db database file

  threshold : float
    edge weight confidence threshold expressed as a value in [0,1]
    include all edges with a confidence greater than or equal to the threshold

  Returns
  -------
  G : nx.Graph
    protein-protein interaction network with Ensembl protein ids as node ids
  """
  # skip header
  line_no = 1
  fh.next()

  G = nx.Graph()
  for line in fh:
    line_no += 1
    p1, p2, conf = parse_line_abc(line)
    if(conf >= threshold):
      G.add_edge(p1, p2, {'weight': conf})
  return G

def string_to_abc(ifh, ofh):
  """
  Convert string database to "abc"-format

  Parameters
  ----------
  ifh : file-like
    String database

  ofh : file-like
    Output file handle to write abc graph to
  """
  # skip header
  ifh.next()
  for line in ifh:
    p1, p2, conf = parse_line_abc(line)
    ofh.write('\t'.join(map(str, [p1, p2, conf])) + '\n')
