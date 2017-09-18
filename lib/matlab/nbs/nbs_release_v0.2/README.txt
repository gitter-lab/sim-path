http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26
http://chianti.ucsd.edu/~mhofree/NBS/nbs_release_v0.2.tgz

Network-based stratification (NBS)
=================================
Version 0.2.00

NBS is a clustering approach designed to uncover cluster from somatic 
mutation data in cancer. It is described in detail in the following 
Publication:

Hofree, M., Shen, J. P., Carter, H., Gross, A., & Ideker, T. (2013). 
Network-based stratification of tumor mutations. Nat Methods. doi: 10.1038/nmeth.2651


This package makes use of several external packages we provide as a curtsey 
in the ./external directory. 

This code was last tested using:
Matlab version: 8.0.0.783 (R2012b)

R:  R version 3.0.0 (2013-04-03)
With the following required packages:
survival
R.matlab
impute
pamr
mclust

System requirements:
===================

Running this package will require a machine with 8 GB of RAM. Consensus 
Clustering is a compute intensive procedure. Much of the analysis we mention in
the above manuscript was performed in parallel on a PBS Torque compute cluster. 
We will release an updated version which will include our framework for 
running the NBS method on a compute cluster in the near future. Please check for 
updates or email us directly.


Installation:
============
1. Unpack project and data files. 
2. Start matlab
2. Mex compile mtimesx for your local environment by running:
cd ./external/mtimesx/
mtimesx_build.m
3. Add library and data directories to your matlab path by changing the 
appropriate variables in the demo files:
library_path
basedata_path


Usage examples:
==============

demo_NBS_simulated.m - Demonstrates simple usage of the NBS method on simulated 
  data.

demo_NBS_TCGA.m - Demonstrates simple usage of the NBS method on publically 
  available uterine data.

NBS_subtypes_to_expression.R - Demonstrates how to use standard pamr (nearest  
shrunken centroid) method in order to translate NBS ovarian subtypes to 
expression subtypes.




===========================================================================
This software is Copyright © 2013 
The Regents of the University of California. All Rights Reserved.

Permission to copy, modify, and distribute this software and its 
documentation for educational, research and non-profit purposes, without 
fee, and without a written agreement is hereby granted, provided that the
 above copyright notice, this paragraph and the following three paragraphs 
appear in all copies.

Permission to make commercial use of this software may be obtained by 
contacting:
Technology Transfer Office
9500 Gilman Drive, Mail Code 0910
University of California
La Jolla, CA 92093-0910
(858) 534-5815
invent@ucsd.edu

This software program and documentation are copyrighted by The Regents of 
the University of California. The software program and documentation are 
supplied "as is", without any accompanying services from The Regents. 
The Regents does not warrant that the operation of the program will be 
uninterrupted or error-free. The end-user understands that the program was 
developed for research purposes and is advised not to rely exclusively on 
the program for any reason.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO
ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING
OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, 
AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO
PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS.
