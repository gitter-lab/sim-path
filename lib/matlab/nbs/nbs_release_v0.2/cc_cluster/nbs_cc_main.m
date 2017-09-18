function [cc_cluster,H_list,sample_data,outlierVect] = ...    
    nbs_cc_main(gind_mat,gene_id_gmap,znetwork,run_options,glapD)
%
% NBS_cc_main - main routine for running consensous clusterling NMF 
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
%   outlierVect - N/A
%
% Matan Hofree (mhofre@cs.ucsd.edu)
% V. 0.2.0 (2/15/2013)
%
% This software is Copyright © 2013 
% The Regents of the University of California. All Rights Reserved.
%
% Please see license.txt for full details.


    % Default settings 
    run_default.nsample = 100;
    run_default.smp_ind = 1;
    run_default.smp_feat = 0.80;    
    run_default.K = 4;         
    run_default.min_indiv = 10; 
    run_default.min_mutations = 9;  
    run_default.propV = 0.6;
    run_default.nmf_type = 'nmf';
    run_default.dis = true;    
    run_default.pb = true;    
    run_default.normalize_rows = 1; 
    run_default.isNetwork = 1;
    
    run_default.cc_type = { 'hard' 'soft' 'euclid' 'corr' 'kldiv' };
    
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
        fprintf(1,'NBS Running with the following options:\n');
        disp(run_options);
    end
    
    % 1 Network NMF
    % Number of indiv/features
    [Nind,Dfeat] = size(gind_mat);
    Nsample = round(Nind*run_options.smp_ind);
    Dsample = round(Dfeat*run_options.smp_feat);                 
    
    if (run_options.pb)
        progressbar();
    end
    % Run clustering
    for iC = 1:run_options.nsample
        tic;               
        % Dsample_vect{iC} = randsample(Dfeat,Dsample);
        % Nsample_vect_t = randsample(Nind,Nsample);
        
        % Just to avoid using the oh so lucrative stats toolbox!
        zperm = randperm(Dfeat);
        Dsample_vect{iC} = zperm(1:Dsample);
        zperm = randperm(Nind);
        Nsample_vect_t = zperm(1:Nsample);

        gene_id_sample = gene_id_gmap(Dsample_vect{iC});
        gind_sample = gind_mat(Nsample_vect_t,Dsample_vect{iC});
        
        % Make sure there are enough mutations
        gind_sample_mutation_count = sum(gind_sample,2);
        gind_sample = gind_sample(gind_sample_mutation_count > run_options.min_mutations,:);
        Nsample_vect{iC} = Nsample_vect_t(gind_sample_mutation_count > run_options.min_mutations);
        
        % Run clustering
        [Tnet{iC},H{iC},outlierVect{iC}] = cluster_data(gind_sample,gene_id_sample,znetwork,run_options,glapD);
        
           
        if (run_options.pb)
            progressbar(iC/run_options.nsample);
        end
        if (run_options.dis)
            fprintf(1,'Round %d time - %f\n',iC,toc);
        end
    end

    fprintf(1,'Done clustering compiling results\n');
    
    
    if (strcmpi(run_options.nmf_type,'nmfOD') ~= 1)
        outlierVect = [];        
    end
    cc_cluster = hmat_to_cc_cluster(run_options.cc_type,run_options.K,Tnet,H,Nsample_vect,Nind,outlierVect);
    H_list = H;
    
    sample_data.Dsample_vect = Dsample_vect;
    sample_data.Nsample_vect = Nsample_vect;
        
    
end 

