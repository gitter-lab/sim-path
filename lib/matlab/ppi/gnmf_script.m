function [ ] = gnmf_script(mat_path, out_w, out_h, out_res, k, wu_init)
% gnmf_script - main method for running graph-regularized nmf
% Input: 
%   mat_path - filepath to a .mat file containing the variables
%     'smoothed_mat' containing the smoothed n_genes x n_screens matrix
%     'laplacian' containing a sparse representation of the laplacian
%   out_w - filepath to write W as a csv file
%   out_h - filepath to write H as a csv file
%   out_res - filepath to write residual of reconstruction error of X with X' := W H
%   k - number of clusters for nmf
%   wu_init - true/false flag to initizalize H with a sample from rows of X
if(nargin < 5)
  error('gnmf_script: mat_path, out_w, out_h, out_res, k are required arguments')
end
if(nargin < 6)
  wu_init = false;
end

% Force mcc to include mtimesx_sparse.m
arr1 = sparse(eye(4));
arr2 = sparse(eye(4) * 2);
rv = mtimesx_sparse(arr1, 'N', arr2, 'N');

%% Process Inputs
if(ischar(k))
  k = str2num(k);
end
load(mat_path);
% TODO verify existence of variables

K = double(laplacian); % MTIMESX cannot multiple int * double
K(~isfinite(K)) = 0; % TODO where did NaNs come from?
% NOTE falpha, fbeta, and fgamma are all deprecated

% nbs requires that n_rows of smoothed_mat are the same as n_rows K
sm_size = size(smoothed_mat);
K_size = size(K);
if(sm_size(2) == K_size(1))
  % then transpose smoothed_mat
  smoothed_mat = smoothed_mat';
end

option.iter = 200;

if(wu_init)
  % provide random initializations for nbs following Wu 2016:
  % sample k rows from smoothed_mat
  [r,c] = size(smoothed_mat);
  rng('shuffle');
  row_inds = unidrnd(r, 1, k);
  row_inds
  Hext = zeros(k, c);
  for i=1:k
    Hext(i,:) = smoothed_mat(row_inds(i), :);
  end
  [W, H, tstats] = nbs_nmf_cluster_network_nmf(smoothed_mat, k, 0, 0, 0, K, option, Hext);
else
  [W, H, tstats] = nbs_nmf_cluster_network_nmf(smoothed_mat, k, 0, 0, 0, K, option);
end

csvwrite(out_w, W);
csvwrite(out_h, H);
csvwrite(out_res, tstats.finalResidual);
end
