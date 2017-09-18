%% NBS demo of recovering and simulating data
% 

library_path = '/cellar/users/mhofree/projects/cancer_strat/code/nbs_release_v0.2';
basedata_path = '/cellar/users/mhofree/projects/cancer_strat/data/nbs_release_v0.2/data';
addpath(genpath(library_path))

%% Load a simulation network
   
baseNetwork = load( [ basedata_path '/networks/PC_NCI_adj_mat.mat'] );
baseClust = load([ basedata_path '/networks/network_clust_PC_NCI.mat'] );

%% Load an NBS propagation network

NBSnetwork = load([ basedata_path '/networks/ST90Q_adj_mat.mat' ]);   

%% Load an NBS graph influence measure 

glap = load([ basedata_path '/networks/glap_subnetwork_ST90.mat' ]);
outDeg = 11;
[knnGlap] = sim_to_knn_glap(glap.glap,outDeg);
clear glap;

%% Load somatic mutation cancer data

baseSMData = load([ basedata_path '/TCGA_somatic_mutations/somatic_data_OV.mat' ]);

min_mutation = 10;
baseSMData.sample_id(sum(baseSMData.gene_indiv_mat,2) < min_mutation) =  [];
baseSMData.gene_indiv_mat((sum(baseSMData.gene_indiv_mat,2) < min_mutation),:) =  [];

%% Simulation clustering to NBSnetwork order

oldidx2id_map = containers.Map(1:length(baseNetwork.network.key), nanvalues(baseNetwork.network.KeytoEntrezIDMap,baseNetwork.network.key));
    
for i = 1:length(baseClust.hq_clus)
    ctval = nanvalues(oldidx2id_map,baseClust.hq_clus{i});
    hq_clus_sub{i} = unique(ctval(~isnan(ctval)));
end
    
%% Run a network based simulation based on NCI cancer pathways

% Simulation parameters:
% Number of subtypes to simualte
optionsDefault.kcluster = 4;
% Number of patients in simulated dataset
optionsDefault.sampleSize =  200;
% Number of genes in simulated dataset
optionsDefault.geneNumber =  8000;
% Number of pathways per cluster to aim for
optionsDefault.numPathWayPerClust = 2;
% Minimum total size of 'cancer pathway'
optionsDefault.pathway_min = 30;
% Maximum total size of 'cancer pathway'
optionsDefault.pathway_max = 50;
% Maximum number of cancer genes reassigned to pathway
% (Maximum number of 'driver' genes)
optionsDefault.p_max = 16;
% Minimum number of cancer genes reassigned to pathway
% (Minimum number of 'driver' genes)
optionsDefault.p_min = 1;
% Proportion of a sample's mutations re-assigned to pathwya
optionsDefault.p_in_pathway = 0.05;
% Number of sub-pathways shared among multiple subtypes
optionsDefault.overlap_clust = 0;
% Minimum size of module to be included in cancer pathway
optionsDefault.minPWsize = 5;

[sim_matrix,sim_indiv_clabel,bkg_gene_ids,bkg_matrix,signal_pathways,kcnt] = ...
    simulate_simple_network_data(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,[],hq_clus_sub,optionsDefault);

sim_matrix =  double(sim_matrix>0);

%% Standard NMF

cczoptions.iter = 1000;
cczoptions.tof = 1e-3;
cczoptions.dis = true;
cczoptions.dist = 'nnls';
cnum = 4;

[A,Y]=nmfnnls(sim_matrix',cnum,cczoptions);

indClust=NMFCluster(Y);
[table, chi2, p] = crosstab(indClust,sim_indiv_clabel);
[AR,RI,MI,HI]=RandIndex(indClust,sim_indiv_clabel);
fprintf('Sim_matrix: P:%e AR:%f\n',p,AR);


%% Run NBS w/o consensous clustering

[indClust] = ...    
    nbs_single_main(sim_matrix,bkg_gene_ids,NBSnetwork.network,[],knnGlap);
[table, chi2, p] = crosstab(indClust,sim_indiv_clabel);
[AR,RI,MI,HI]=RandIndex(indClust,sim_indiv_clabel);
fprintf('Expgeno: P:%e AR:%f\n',p,AR);


%% Run consensous clustering (CC) NBS

% Type of NMF to use:
propnmf_options.nmf_type = 'netnmf';
% Number of times to perform NMF as part of CC
% (100 is probably the bare minimum)
propnmf_options.nsample = 100;
% Proportion of rows (patient samples) to include in CC
propnmf_options.smp_ind = 0.8;
% Proportion of columns (patient samples) to include in CC
propnmf_options.smp_feat = 0.8;
% Number of clusters 
propnmf_options.K = [ 4 ];
% Minimum number of mutations persample to be includeded 
propnmf_options.min_mutations = 9;
% Propagation value
propnmf_options.proV = NBSnetwork.network.propVal;
% Verbose output
propnmf_options.dis = true;
% Show progress bar (requries X11)
propnmf_options.pb = true;
% Quantile normzlize
propnmf_options.normalize_rows = 1;
% Use netwpork
propnmf_options.isNetwork = 1;

[NBS_cc_cluster] = ...    
    nbs_cc_main(sim_matrix,bkg_gene_ids,NBSnetwork.network,propnmf_options,knnGlap);

[NBS_cc_label] = calc_label(NBS_cc_cluster.clust_hard.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(NBS_cc_label,sim_indiv_clabel);
[AR,RI,MI,HI]=RandIndex(NBS_cc_label(~isnan(NBS_cc_label)),sim_indiv_clabel(~isnan(NBS_cc_label)));
fprintf('Expgeno: P:%e AR:%f\n',p,AR);


%% Run standard consensous clustering NMF

propnmf_options.nmf_type = 'nmf';
propnmf_options.isNetwork = 0;

[ctrl_cc_cluster] = ...    
    nbs_cc_main(sim_matrix,bkg_gene_ids,[],propnmf_options);

[ctrl_cc_label] = calc_label(ctrl_cc_cluster.clust_soft.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(ctrl_cc_label,sim_indiv_clabel);
[AR,RI,MI,HI]=RandIndex(ctrl_cc_label(~isnan(ctrl_cc_label)),sim_indiv_clabel(~isnan(ctrl_cc_label)));
fprintf('Expgeno: P:%e AR:%f\n',p,AR);


