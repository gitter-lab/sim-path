function [t_cluster,H_list,outlierVect] = ...    
    nbs_single_main(gind_mat,gene_id_gmap,znetwork,run_options,knnGlap)
%
% NBS_single_main - main routine for running consensous clusterling NMF 
% Input:
%   gind_mat - A binary patinet x somatic mutation matrix 
%   gene_id_gmap - A vector of entrez IDs identifying the mutation in 
%        gind_mat matrix
%   znetwork - A preprocessed network in the form of a propNetwork object 
%   run_options - Input options
%   glapD - A distance functions on the network
%
% Output:
%   cc_cluster - A cell array of out consensous clusering arrays
%   H_list - Output soft assignment matrix 
%   sample_Data - Samples used in each cc sample
%   outlierVect - (outlier detecting NMF only)
%
% Matan Hofree (mhofre@cs.ucsd.edu)
% V. 0.2.0 (2/15/2013)
%
% This software is Copyright © 2013 
% The Regents of the University of California. All Rights Reserved.
%
% Please see license.txt for 


% Default settings
run_default.K = 4;
run_default.min_indiv = 10;
run_default.min_mutations = 9;
run_default.propV = 0.7;
run_default.nmf_type = 'netnmf';
run_default.dis = true;
run_default.pb = true;
run_default.normalize_rows = 1;
run_default.isNetwork = 1;
    
% Defaule clustering settings
run_default.zoptions.iter = 1500;
run_default.zoptions.gamma = 200;
run_default.zoptions.tof = 1e-4;
run_default.zoptions.dis = false;
run_default.zoptions.distance = 'nnls';
    
if (exist('run_options','var') == 1)
    run_options = mergeOption(run_options,run_default);
else
    run_options = run_default;
end

if (exist('glapD','var') == 0)
    glapD = [];
end

if (run_options.dis)
    fprintf(1,'Single round NBS running with the following settings:\n');
    disp(run_options);
end
  
[t_cluster,H_list,outlierVect] = cluster_data(gind_mat,gene_id_gmap,znetwork,run_options,knnGlap);
  
