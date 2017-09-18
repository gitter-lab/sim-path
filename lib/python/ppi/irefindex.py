import sys

import irefindex_parser
reload(irefindex_parser)
from irefindex_parser import *

import metrics_nx
reload(metrics_nx)
from metrics_nx import *

try:
  import metrics_gt
  reload(metrics_gt)
except ImportError:
  sys.stderr.write("[warning] Cannot import graph_tool\n")
