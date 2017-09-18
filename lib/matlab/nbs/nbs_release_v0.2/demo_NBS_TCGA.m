%% NBS demo of recovering and simulating data
% 

library_path = '/cellar/users/mhofree/projects/cancer_strat/code/nbs_release_v0.2';
basedata_path = '/cellar/users/mhofree/projects/cancer_strat/data/nbs_release_v0.2/data';
addpath(genpath(library_path))

%% Load an NBS propagation network

NBSnetwork = load([ basedata_path '/networks/ST90Q_adj_mat.mat' ]);   

%% Load an NBS graph influence measure 

glap = load([ basedata_path '/networks/glap_subnetwork_ST90.mat' ]);
outDeg = 11;
[knnGlap] = sim_to_knn_glap(glap.glap,outDeg);
clear glap;

%% Load somatic mutation cancer data

baseSMData = load([ basedata_path '/TCGA_somatic_mutations/somatic_data_UCEC.mat' ]);

min_mutation = 10;
baseSMData.sample_id(sum(baseSMData.gene_indiv_mat,2) < min_mutation) =  [];
baseSMData.gene_indiv_mat((sum(baseSMData.gene_indiv_mat,2) < min_mutation),:) =  [];

baseSMData.sample_id = regexprep(baseSMData.sample_id,'(TCGA-..-....).*','$1');

%% Load phenotype

load([basedata_path '/TCGA_somatic_mutations/UCEC_clinical_phenotype.mat' ]);

%% 

[~,ia,ib] = intersect(baseSMData.sample_id,UCECppheno.sample_id);

jointPheno = cellfun(@(x,y)([x '-' y]),UCECppheno.hist_simple,UCECppheno.grade_simple,'uniformoutput',0);
cphenoJoint.sample_id = UCECppheno.sample_id(ib);
cphenoJoint.phenoMat = [grp2idx(jointPheno(ib)) grp2idx(UCECppheno.hist_simple(ib)) grp2idx(UCECppheno.grade_simple(ib)) ]


%% Standard NMF

cczoptions.iter = 1000;
cczoptions.tof = 1e-3;
cczoptions.dis = true;
cczoptions.dist = 'nnls';
cnum = 3;

[A,Y]=nmfnnls(baseSMData.gene_indiv_mat',cnum,cczoptions);

indClust_ctrl=NMFCluster(Y);
[table, chi2, p] = crosstab(indClust_ctrl(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(indClust_ctrl(ia),cphenoJoint.phenoMat(:,1));
fprintf('Sim_matrix: P:%e AR:%f\n',p,AR);


%% Run NBS w/o consensous clustering

% Type of NMF to use:
propnmf_options.nmf_type = 'netnmf';
% Number of clusters 
propnmf_options.K = [ 3 ];
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

[indClust] = ...    
    nbs_single_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,NBSnetwork.network,propnmf_options,knnGlap);
%%
[table, chi2, p] = crosstab(indClust(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(indClust(ia),cphenoJoint.phenoMat(:,1));
fprintf('NBS: P:%e AR:%f\n',p,AR);
    


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
propnmf_options.K = [ 3 4 ];
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
    nbs_cc_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,NBSnetwork.network,propnmf_options,knnGlap);
%%
[NBS_cc_label] = calc_label(NBS_cc_cluster.clust_hard.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(NBS_cc_label(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(NBS_cc_label(ia),cphenoJoint.phenoMat(:,1));
fprintf('CC NBS: P:%e AR:%f\n',p,AR);


%% Run standard consensous clustering NMF

propnmf_options.nmf_type = 'nmf';
propnmf_options.isNetwork = 0;

[ctrl_cc_cluster] = ...    
    nbs_cc_main(baseSMData.gene_indiv_mat,baseSMData.gene_id_all,[],propnmf_options);

[ctrl_cc_label] = calc_label(ctrl_cc_cluster.clust_hard.cc_matrix{cnum},cnum,'av',1);
[table, chi2, p] = crosstab(ctrl_cc_label(ia),cphenoJoint.phenoMat(:,1));
[AR,RI,MI,HI]=RandIndex(ctrl_cc_label(ia),cphenoJoint.phenoMat(:,1));
fprintf('Standard CC NMF: P:%e AR:%f\n',p,AR);

