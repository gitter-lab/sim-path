"""
Functions for working with EBI's HGNC database
"""
# origin 0 field indexes
HGNC_ID = 0 # HGNC:5
HGNC_SYMBOL = 1 # A1BG
DATE_APPROVED_RESERVED = 14 # YYYY-MM-DD
DATE_SYMBOL_CHANGED = 15 # YYYY-MM-DD

def parse_hgnc_db(db_fh, n_lines=None):
  rv = {}
  # advance past header
  db_fh.next()
  if n_lines is None:
    n_lines = float("inf")
  line_no = 0
  for line in db_fh:
    line = line.rstrip()
    tokens = line.split("\t")
    id_str = tokens[HGNC_ID]
    hgnc_label, hgnc_id_str = id_str.split(":")
    hgnc_symbol = tokens[HGNC_SYMBOL]
    rv[hgnc_id_str] = hgnc_symbol
    line_no += 1
    if(line_no > n_lines):
      break
  return rv
