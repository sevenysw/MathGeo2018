function [filter_bank] = filter_learning_3D(data, lambda, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning orthogonal filter bank  from input
% input:
%   data		-	data is a matrix with p predictors and n obervations
%   lambda		-	parameter for thresholding
%   (optional)
%   opts:
%		A		-	initialization for A, whose atoms are orthogonal with each other.
%   	nIter	-	number of Iteration for learning, default value (=50)
% 
% output:
%   filter_bank	-	output of data-driven filter bank
%
%Reference: Jian-feng Cai, H. Ji, Z. Shen and Guibo Ye,  
%Data-driven tight frame construction and image denoising ,
%Applied and Computational Harmonic Analysis, 37 (1), 89-105, Jul. 2014
%
%Author: Chenlong Bao, Yuhui Quan, Sibin Huang, and Hui Ji
%
%Last Revision: 25-May-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Exception Management

if nargin < 2; error('No enough input!'); end
%% default options
if ~exist('opts','var')   	opts = struct;   	   end
if isfield(opts,'A')     	A = opts.A;      	   else  A = [];       end
if isfield(opts,'nIter')    nIter = opts.nIter;    else  nIter = 50;   end

%% whether A is empty 
flag_A = isempty(A);

%% Input and initialization 
data = double(data);   	  % input data
p    = size(data, 1);     % dimension of data
if ~flag_A
	A = double(A);		  
	tmpMat = (eye(p) - A*A')*data;
end
% D = DctInit(round(p^(1/2))); 
[D,~] = gen_dic_by_iwt_3d(round(p^(1/3)),'haar');
if ~flag_A
	r = size(A, 2);
	rperm = randperm(p);
	D = D(:,rperm(1:p-r));
end
D = double(D);			 % input dictionary D

%% Learning loop
iter = 0;
while iter < nIter
    % Soft thresholding for solving L-0 minimization
	tmpCoef = wthresh(D'*data, 'h', lambda);
    % SVD for solving dictionary update
	if ~flag_A
		[U, S, V] = svd(tmpMat*tmpCoef',0);
    else

		[U, S, V] = svd(data*tmpCoef',0);
    end
    D = U*V';
    % next loop
    iter = iter + 1;
end

%% Output
filter_bank = [A, D];      %Dictionary