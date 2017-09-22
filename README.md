# sim-path

sim-path is a software project to improve subnetwork identification with pathway simulation. 
To use the software, run the script @pipeline\_simulation.py@ with the required arguments.
For example,

```bash
./pipeline_simulation.py ~/data/irefindex/9606.mitab.04072015.txt ~/data/pathways/ ~/data/simpath_out/
```

The required arguments include
1. A background protein-protein interaction network
2. A set of pathways from which to simulate genetic screen hits
3. An output directory to place results

It is our intent to support different protein-protein interaction databases.
[iRefIndex](http://irefindex.org/download/irefindex/data/archive/release_14.0/psi_mitab/MITAB2.6/9606.mitab.07042015.txt.zip)
is the recommended protein-protein interaction database to use.
It is possible to also use the [STRING](https://string-db.org/download/protein.links.full.v10.5/9606.protein.links.full.v10.5.txt.gz)
protein-protein association database.
However, doing so requires additional arguments to map identifiers and is not yet documented.

Our experiments include reference pathways from the [Kyoto Encyclopedia of Genes and Genomes (KEGG)](http://www.genome.jp/kegg/).
We provide a high-level interface that connects to the KEGG REST API to download reference pathways in @kegg.R@.
For our experiments, we have downloaded all pathways listed by @kegg.R -l@
We chose to use KEGG pathways in our simulations because they are more representative of real biological processes than randomly generated subnetworks of a global protein interaction network.
However, any set of biological pathways may be used with this software.
These pathways must be represented in the @.graphml@ file format and contain node identifiers that match those in the background protein-protein interaction network.
Note that the KEGG REST API for downloading pathways is intended for academic use only.

Results will be placed in the output directory provided.

# Requirements

* Python 2.7
  * numpy
  * scikit-learn
  * networkx 
* R 3.4.0
  * KEGGREST
  * KEGGgraph
  * igraph

# License

This repository contains original work by Aaron Baker and redistributed work from the authors 
of the @nbs@ and @stratipy@ projects. Each of these works are licensed separately.
Unless otherwise noted, all files in this repository are authored by Aaron Baker and licensed
according to the MIT license contained in the file LICENSE.txt.
The first of two exceptions to this is the redistributed files of the @nbs@ project located at
@lib/matlab/nbs/nbs\_release\_v0.2/@, which are licensed according to LICENSE-NBS.txt.
The second exception is the redistributed files of the @stratipy@ project located at
@lib/python/stratipy/@, which are licensed according to LICENSE-STRATIPY.txt
