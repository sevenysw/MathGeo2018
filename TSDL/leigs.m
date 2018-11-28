function [E,V] = leigs(DATA, TYPE, PARAM, NE)  

% Laplacian Eigenmaps Algorithm 
%
% please refer to University of Chicago 
% Computer Science Technical Report TR-2002-01
% Mikhail Belkin, Partha Niyogi
% Laplacian Eigenmaps for Dimensionality Reduction and Data Representation
% Note that Laplacian, not normalized Laplacian is used here
% http://www.cs.uchicago.edu/research/publications/techreports/TR-2002-1
%
%
% Calculate the graph laplacian of the adjacency graph of data set DATA.
%
% L = laplacian(DATA, TYPE, PARAM, NE)  
% 
%   DATA - NxK matrix. Data points are rows. 
%   TYPE - string 'nn' or string 'epsballs'
%   PARAM - integer if TYPE='nn', real number if TYPE='epsballs'
%   NE - number of eigenvectors
%
% Returns: 
%   E - NxNE matrix with eigenfunctions, 
%   V is NExNE matrix with eigenvalues on the diagonal
%
% Author: 
%
%   Mikhail Belkin 
%   misha@math.uchicago.edu
%

L = laplacian(DATA, TYPE, PARAM);

% normalized Laplacian
% for i=1:size (L)
%   D(i,i) = sum(L(i,:));
%   if (D(i,i) ~= 0)
%      DD(i,i) = 1/sqrt(D(i,i));
%   else disp ('warning 0');
%      DD(i,i) = 0;
%   end
% end
% LL=DD*(D-L)*DD;


opts.tol = 1e-9;
opts.issym=1; 
opts.disp = 5; 
[E,V] = eigs(L,NE,'sm',opts);

%A = DD*A;

 
