def detect_db(fh):
  """
  Parameters
  ----------
  fh : file-like object

  Returns
  -------
  db_type : str
    one of 'irefindex', 'string'

  Notes
  -----

  STRING
  ======
  head -n2 9606.protein.links.full.v10.5.txt
  protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score
  9606.ENSP00000000233 9606.ENSP00000263431 0 0 0 0 0 0 53 0 176 0 0 0 128 260

  iRefIndex
  ========
  head -n2 9606.mitab.04072015.txt
  #uidA	uidB	altA	altB	aliasA	aliasB	method	author	pmids	taxa	taxb	interactionType	sourcedb	interactionIdentifier	confidence	expansion	biological_role_A	biological_role_B	experimental_role_A	experimental_role_B	interactor_type_A	interactor_type_B	xrefs_A	xrefs_B	xrefs_Interaction	Annotations_A	Annotations_B	Annotations_Interaction	Host_organism_taxid	parameters_Interaction	Creation_date	Update_date	Checksum_A	Checksum_B	Checksum_Interaction	Negative	OriginalReferenceA	OriginalReferenceB	FinalReferenceA	FinalReferenceB	MappingScoreA	MappingScoreB	irogida	irogidb	irigid	crogida	crogidb	crigid	icrogida	icrogidb	icrigid	imex_id	edgetype	numParticipants
  uniprotkb:A0A024R3E3	uniprotkb:A0A024R3E3	entrezgene/locuslink:335|genbank_protein_gi:4557321|refseq:NP_000030|rogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606|irogid:122812673	entrezgene/locuslink:335|genbank_protein_gi:4557321|refseq:NP_000030|rogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606|irogid:122812673	hgnc:APOA1|uniprotkb:A0A024R3E3_HUMAN|uniprotkb:APOA1_HUMAN|crogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606|icrogid:122812673	hgnc:APOA1|uniprotkb:A0A024R3E3_HUMAN|uniprotkb:APOA1_HUMAN|crogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606|icrogid:122812673	-	-	pubmed:9003180|pubmed:9200714|pubmed:9356442	taxid:9606(Homo sapiens)	taxid:9606(Homo sapiens)	-	MI:0462(bind)	bind:75986|rigid:WRUQaMHXGmnC/H/BzyolIfyaa7Y|edgetype:X	hpr:141818|lpr:1|np:8	none	MI:0000(unspecified)	MI:0000(unspecified)	MI:0000(unspecified)	MI:0000(unspecified)	MI:0326(protein)	MI:0326(protein)	-	-	-	-	-	-	-	-	2015-04-07	2015-04-07	rogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606	rogid:yEzPDeU8Uu/43dkLLOBAy6ey1vs9606	rigid:WRUQaMHXGmnC/H/BzyolIfyaa7Y	false	GenBank:NP_000030	GenBank:NP_000030	refseq:NP_000030	refseq:NP_000030	PD	PD	122812673	122812673	881764	yEzPDeU8Uu/43dkLLOBAy6ey1vs9606	yEzPDeU8Uu/43dkLLOBAy6ey1vs9606	WRUQaMHXGmnC/H/BzyolIfyaa7Y	122812673	122812673	881764	-	X	2
  """
  string_header = """protein1 protein2 neighborhood neighborhood_transferred fusion cooccurence homology coexpression coexpression_transferred experiments experiments_transferred database database_transferred textmining textmining_transferred combined_score""".split()
  irefindex_header = """#uidA uidB  altA  altB  aliasA  aliasB  method  author  pmids taxa  taxb  interactionType sourcedb  interactionIdentifier confidence  expansion biological_role_A biological_role_B experimental_role_A experimental_role_B interactor_type_A interactor_type_B xrefs_A xrefs_B xrefs_Interaction Annotations_A Annotations_B Annotations_Interaction Host_organism_taxid parameters_Interaction  Creation_date Update_date Checksum_A  Checksum_B  Checksum_Interaction  Negative  OriginalReferenceA  OriginalReferenceB  FinalReferenceA FinalReferenceB MappingScoreA MappingScoreB irogida irogidb irigid  crogida crogidb crigid  icrogida  icrogidb  icrigid imex_id edgetype  numParticipants""".split()

  rv = None
  fh.seek(0)
  fh_header = fh.next().rstrip().split()
  if fh_header == string_header:
    rv = "string"
  elif fh_header == irefindex_header:
    rv = "irefindex"
  else:
    raise RuntimeError("Unrecognized database file: it is not one of iRefIndex, STRING")
  fh.seek(0)
  return rv
