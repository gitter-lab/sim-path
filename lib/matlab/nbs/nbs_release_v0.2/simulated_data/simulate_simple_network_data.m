function [sim_matrix,sim_indiv_clabel,bkg_gene_ids,bkg_matrix,signal_pathways,kClusterCounter] = ... 
    simulate_simple_network_data_cluster_patient_distrib_minmax(sample_data,sample_id,simnetwork,pathway_set,simoptions)
% Simulate a cancer cohort with signal embedded to a set of network modules
%
% Matan Hofree (mhofre@cs.ucsd.edu)
% V. 0.2.0 (2/15/2013)
%
% This software is Copyright © 2013 
% The Regents of the University of California. All Rights Reserved.
%
% Please see license.txt for full details.


    optionsDefault.kcluster = 4;
    optionsDefault.sampleSize =  400;
    optionsDefault.geneNumber =  5000;
    optionsDefault.numPathWayPerClust = 2;
    optionsDefault.scrambleGeneID = true;
    optionsDefault.pathway_min = 10;
    optionsDefault.pathway_max = 120;
    optionsDefault.p_max = 0.9;
    optionsDefault.p_min = 5;  
    optionsDefault.p_in_pathway = 0.2; 
    optionsDefault.rand_choose = 1;   
    optionsDefault.overlap_clust = 0;
    optionsDefault.minPWsize = 5;
    if nargin<4
        simoptions=optionsDefault;
    else
        simoptions=mergeOption(simoptions,optionsDefault);
    end        

    disp(simoptions);

    sIdx = randperm(length(pathway_set));
    pathway_set = pathway_set(sIdx);    
    
    % Assign pathways to clustersS
    [path_assign,kClusterCounter] = pathway_assign(simoptions,pathway_set);
    
    pathway_set_selected = pathway_set(concatCell(path_assign));    
    pathway_genes = unique(concatCell(pathway_set_selected));
    
    % Generate background data
    [bkg_matrix,bkg_gene_ids] = generate_bkg(sample_data,sample_id,pathway_genes,simoptions);        
    
    % Imbed mutations with a higher freq then background    
    [sim_matrix,bkg_gene_ids,sim_indiv_clabel,signal_pathways] = simulate_bkg_with_pathways_patient(bkg_matrix,bkg_gene_ids,pathway_set,path_assign,simoptions);    

end

function [path_assign,kClusterCounter] = pathway_assign(simoptions,pathway_set)
    
    kClusterCounter = zeros(simoptions.kcluster,1);
    kClusterSum = zeros(simoptions.kcluster,1);   
    pwsize = cellfun(@(x)(length(x)),pathway_set);    
    path_assign = {};
    
    % Assign pathways
    kfactor = simoptions.numPathWayPerClust;
    minSize = max(simoptions.minPWsize,floor(simoptions.pathway_max/kfactor));
    
    
    kNeedChange = (kClusterSum < simoptions.pathway_min | kClusterSum > simoptions.pathway_max);
    
    loopD = 0;
    while (any(kNeedChange))
        loopD = loopD + 1;
        for i = find(kNeedChange)'          
            candidatePathways = find(pwsize > minSize);
            while ( isempty(candidatePathways) || loopD > simoptions.numPathWayPerClust*2 )
                kfactor = kfactor + 1;
                minSize = max(simoptions.minPWsize,floor(simoptions.pathway_max/kfactor));
                candidatePathways = find(pwsize > minSize);
                loopD = 0;
            end            
                
            
            cnt_p = randsample(candidatePathways,1);
            if ( kClusterSum(i) > simoptions.pathway_max ) % Switch a pathway
                zselect = randsample(kClusterCounter(i),1);
                repPos = path_assign{i}(zselect);
                
                pwsize(repPos) = length(pathway_set{repPos});
                kClusterSum(i) = kClusterSum(i) - pwsize(repPos);
                
                path_assign{i}(zselect) = cnt_p;                
                kClusterSum(i) = kClusterSum(i) + pwsize(cnt_p);
                pwsize(cnt_p) = 0;                
            else % Add a pathway
                kClusterCounter(i) = kClusterCounter(i) + 1;
                if ( isempty(path_assign) || length(path_assign) < i)
                     path_assign(i) = { cnt_p };
                else
                    path_assign{i}(kClusterCounter(i)) = cnt_p;
                end
                kClusterSum(i) = kClusterSum(i) + pwsize(cnt_p);
                pwsize(cnt_p) = 0;                
            end
        end
        kNeedChange = (kClusterSum < simoptions.pathway_min | kClusterSum > simoptions.pathway_max);
    end    
        
    % Generate overlap
    overlap_counter = ones(simoptions.kcluster,1)*simoptions.overlap_clust;
    overlap_counter = min(overlap_counter,kClusterCounter);
    
    %isOverlap = false(simoptions.kcluster,simoptions.numPathWayPerClust);    
    isOverlap = cellfun(@(x)false(size(x)),path_assign,'uniformoutput',false);
    
    overlap_counter = min(overlap_counter,kClusterCounter-1);
    
    while any(overlap_counter(i) > 0)
        for i = 1:simoptions.kcluster             
            if (overlap_counter(i) > 0)
                overlap_counter(i) = overlap_counter(i) - 1;                 
                
                % Find a resonable switch 
                % while (start || (isOverlap(clustI,clustJ) || isOverlap(i,clustJ)) )
                clustExI = randsample(setdiff(1:simoptions.kcluster,i),1);
                clustExJ = randsample(length(path_assign{clustExI}),1);
                if ( all(~isOverlap{i} ) ) % Try find empty
                    clustJ = randsample(find(~isOverlap{i}),1);
                else                       % Or choose one at random                    
                    clustJ = randsample(length(path_assign{i}),1);
                end
                
                isOverlap{clustExI}(clustExJ) = 1;
                isOverlap{i}(clustJ) = 1;
                                
                path_assign{i}(clustJ) = path_assign{clustExI}(clustExJ);    
            end
        end
    end

end

function [sim_matrix,bkg_gene_ids,sim_indiv_clabel,p_matrix] = simulate_bkg_with_pathways_patient(bkg_matrix,bkg_gene_ids,pathway_set_good,path_assign,simoptions)
        
    [Nsample_data,Dsample_data] = size(bkg_matrix);
    sim_matrix = bkg_matrix;    
    
    bkg2idx_map = containers.Map(bkg_gene_ids,1:length(bkg_gene_ids));        
    
    % Assign to clusters
    cluster_size = floor(Nsample_data/simoptions.kcluster);    
    sim_indiv_clabel = nan(Nsample_data,1);
    
    for i = 1:simoptions.kcluster
        sim_indiv_clabel(((i-1)*cluster_size+1):(i*cluster_size)) = i;
        
        % p_set = pathway_set_good(path_assign(i,:));
        p_set_all{i} = concatCell(pathway_set_good(path_assign{i}));
        p_matrix{i} = nanvalues(bkg2idx_map,p_set_all{i});
    end
    sim_indiv_clabel((i*cluster_size):end) = i;            
        
    % Mutate with higher freq within pathway
    bkg_matrix_freq = sum(bkg_matrix)/size(bkg_matrix,1);
    
    for i = 1:Nsample_data               
        ps_set_mapped = p_matrix{sim_indiv_clabel(i)};
        
        c_mutation_set = find(sim_matrix(i,:));
        c_cnt_mutations = length(c_mutation_set);
        c_cnt_pathway_mutations = max(simoptions.p_min,ceil(length(c_mutation_set)*simoptions.p_in_pathway));
        
        if (simoptions.p_max > 1)
            c_cnt_pathway_mutations = min(c_cnt_pathway_mutations,simoptions.p_max);            
        end
        
        % Remove original mutations
        reassign_mutation = randperm(c_cnt_mutations,min(c_cnt_pathway_mutations,c_cnt_mutations));
        sim_matrix(i,c_mutation_set(reassign_mutation)) = 0;
        
        % Re-assign mutations by pathway
        pathway_mutations = randperm(length(ps_set_mapped),min(c_cnt_pathway_mutations,length(ps_set_mapped)));
        sim_matrix(i,ps_set_mapped(pathway_mutations)) = 1;                       
    end
end


function [bkg_sim_mat,bkg_sample_gene_id] = generate_bkg(sample_data,sample_gene_id,pathway_genes,simoptions)

    [Nsample_data,Dsample_data] = size(sample_data);
    assert(Dsample_data >= simoptions.geneNumber,'More genes than input');
    pwgenesD = length(pathway_genes);
    
    % Choose a random sample of genes for simulation             
    if (simoptions.scrambleGeneID) 
        fprintf(1,'Note: using randomly assigned gene distributions\n');  
        scr = randsample(Dsample_data,simoptions.geneNumber,true);
        gene_id_diff = setdiff(sample_gene_id,pathway_genes);
        
        scr_label = randperm(length(gene_id_diff));
        
        bkg_sample_gene_id = [gene_id_diff(scr_label(1:(simoptions.geneNumber-pwgenesD))); pathway_genes]';
        bkg_sample_data = sample_data(:,scr);
    else % Use the actual gene parameters for simulation
       fprintf(1,'Note: using original gene distributions\n');
       [~,gene_id_diff_idx] = setdiff(sample_gene_id,pathway_genes);
       [~,gene_id_int] = intersect(sample_gene_id,pathway_genes);
       fprintf(1,'Note: %d of %d pathwway genes found\n',length(gene_id_int),pwgenesD);              
       scr_label = randperm(length(gene_id_diff_idx));
       
       src_idx = [gene_id_diff_idx(scr_label(1:(simoptions.geneNumber-length(gene_id_int)))) gene_id_int];
       
       bkg_sample_gene_id = sample_gene_id(src_idx);
       bkg_sample_data = sample_data(:,src_idx);
    end
        
    bkg_sim_mat = zeros(simoptions.sampleSize,simoptions.geneNumber);
    [~,Dsim_data] = size(bkg_sim_mat);
    
    % Sample genes according to original patient  distributions 
    for i = 1:simoptions.sampleSize
         g_picks = randi(Nsample_data,1);         
         zperm = randperm(Dsim_data);
         bkg_sim_mat(i,:) = bkg_sample_data(g_picks,zperm);
    end
     
    % Sanity check: Test for diviation from original freq.
    bkg_sample_data_freq = sum(bkg_sample_data)/size(bkg_sample_data,1);
    bkg_sim_freq = sum(bkg_sim_mat)/size(bkg_sim_mat,1);
    
    fprintf(1,'Deviation of sample from org data %f\n',norm(bkg_sample_data_freq-bkg_sim_freq))
    
end

