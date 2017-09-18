import graph_tool as gt
import graph_tool.search
from graph_tool import GraphView
from .irefindex import *

class IrefGraph():
  """Store data associated with the iRefIndex database file in this object
  to provide a more convenient interface to network analysis functions.
  """

  def __init__(self, fh, n_rows=None):
    G, rogid_to_vertex = get_graph_tool_ppi(fh, n_rows=None, hgnc=True)
    self.G = G
    self.rogid_to_vertex = rogid_to_vertex
    self.alias_map = get_aliases(fh)

  def enum_deep_neighbors(self, hgncs, depth=1):
    return enum_deep_neighbors(self.G, self.rogid_to_vertex, hgncs, self.alias_map, depth)

  def induce_deep_neighbor_subgraph(self, hgncs, depth=1):
    return induce_deep_neighbor_subgraph(self.G, self.rogid_to_vertex, self.alias_map, hgncs, depth)

def get_graph_tool_ppi(fh, n_rows=None, hgnc=False):
  """Parse graph_tool ppi network from irefindex

  Returns
  -------
  G : graph_tool.Graph
    The parsed irefindex PPI network with a vertex for each unique rogid

  rogid_to_vertex : dict
    mapping of rogid to vertex object in G

  TODO
  ----
  There are many cases where there are multiple HGNCs for a single interactor
  This step takes about 2.5 minutes -- can this be improved by serialize/deserialize in graph tool native format
  """

  # TODO currently too many edges
  try:
    fh.seek(0)
  except AttributeError:
    pass

  header = parse_header(fh)

  if n_rows is None:
    # then parse all data rows
    n_rows = np.Inf

  # Construct interaction graph
  G = gt.Graph(directed=False)
  vprop_rogid = G.new_vertex_property("string")
  if(hgnc):
    vprop_hgnc = G.new_vertex_property("string")
  # TODO shouldnt the property map do this??
  rogid_to_vertex = {}
  edges = set()
  i = 0
  for line in fh:
    line = line.strip()
    datum = line.split("\t")
    if (not is_human_human(header, datum)):
      # include only human-human interactions
      continue
    id_a, id_b = get_checksums(header, datum)
    id_a = id_a[1]
    id_b = id_b[1]

    # create vertices and associate them with properties
    v_a = None
    if(id_a not in rogid_to_vertex):
      v_a = G.add_vertex()
      vprop_rogid[v_a] = id_a[1]
      rogid_to_vertex[id_a] = v_a
    else:
      v_a = rogid_to_vertex[id_a]
    v_b = None
    if(id_b not in rogid_to_vertex):
      v_b = G.add_vertex()
      vprop_rogid[v_b] = id_b[1]
      rogid_to_vertex[id_b] = v_b
    else:
      v_b = rogid_to_vertex[id_b]

    # Graph_Tool doesnt appear to have a different way to have a non multigraph
    # so we do so here
    edge = (v_a, v_b)
    edge_r = (v_b, v_a)
    if not (edge in edges or edge_r in edges):
      G.add_edge(v_a, v_b)
      edges.add(edge)
      edges.add(edge_r)

    if (hgnc):
      # then add hgnc as a node attribute
      # TODO may there be multiple hgnc identifiers?
      hgncs_a, hgncs_b = get_hgnc_pair(header, datum)
      attrs_a = { HGNCS_KEY: hgncs_a }
      attrs_b = { HGNCS_KEY: hgncs_b }
      if(len(hgncs_a) > 0):
        vprop_hgnc[v_a] = hgncs_a[0]
      if(len(hgncs_b) > 0):
        vprop_hgnc[v_b] = hgncs_b[0]

    i += 1
    if i > n_rows:
      break

  # internalize the property maps
  G.vertex_properties['rogid'] = vprop_rogid
  if(hgnc):
    G.vertex_properties['hgnc'] = vprop_hgnc

  return G, rogid_to_vertex

def enum_deep_neighbors(G, rogid_to_vertex, alias_map, hgncs, depth=1):
  """Enumerate the neighborhood of each vertex whose HGNC identifier is given in <hgncs>

  The classical neighborhood N_G(v) = {v' in V | (v, v') in E} is associated with depth=1
  Define a deep neighborhood of depth 2 by N_G(N_G(v)) or N_G^n(v) for n=2

  Returns
  -------
  rv : dict
    mapping of vertex to its deep neighborhood, a list of vertex objects
  """
  rv = {}

  hgnc_alias_map = rekey_aliases(alias_map, 'hgnc')
  for hgnc in hgncs:
    if hgnc not in hgnc_alias_map:
      sys.stderr.write("Unknown hgnc {}\n".format(hgnc))
      continue
    rogid = hgnc_alias_map[hgnc]['rogid']
    vertex = rogid_to_vertex[rogid]

    # new property map and visitor for each BFS
    deep_neighbors = []
    dist_map = G.new_vertex_property("int")
    edge_iter = graph_tool.search.bfs_iterator(G, source=vertex)
    for e in edge_iter:
      dist = dist_map[e.source()] + 1
      if(dist > depth):
        break
      deep_neighbors.append(e.target())
      dist_map[e.target()] = dist
    rv[vertex] = deep_neighbors

  return rv

def induce_subgraph(G, rogid_to_vertex, alias_map, hgncs):
  """Return subgraph induced by vertices defined by hgncs
  """
  vfilt = G.new_vertex_property('bool')
  hgnc_alias_map = rekey_aliases(alias_map, 'hgnc')
  for hgnc in hgncs:
    if hgnc not in hgnc_alias_map:
      continue
    rogid = hgnc_alias_map[hgnc]['rogid']
    vertex = rogid_to_vertex[rogid]
    vfilt[vertex] = True
  sub_G = GraphView(G, vfilt)
  return sub_G

def induce_deep_neighbor_subgraph(G, rogid_to_vertex, alias_map, hgncs, depth=1):
  if 'hgnc' not in G.vertex_properties:
    raise Exception('Missing vertex property \"hgnc\"')

  v_to_ns = enum_deep_neighbors(G, rogid_to_vertex, alias_map, hgncs, depth)
  vs = set()
  for v, ns in v_to_ns.iteritems():
    vs.add(v)
    for n in ns:
      vs.add(n)

  hgncs_neigh = []
  for v in vs:
    hgnc = G.vertex_properties['hgnc'][v]
    if len(hgnc) != 0:
      hgncs_neigh.append(hgnc)

  G_sub = induce_subgraph(G, rogid_to_vertex, alias_map, hgncs_neigh)
  return G_sub
