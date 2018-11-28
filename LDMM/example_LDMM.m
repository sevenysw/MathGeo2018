% Noise attenuation in a low dimensional manifold
% Author: Siwei Yu
% Reference: Yu S, et. al.  Noise attenuation in a low dimensional
% manifold, accepted by geophysics.
% Email:  jma@hit.edu.cn
% March 22, 2017

clear;
load f.mat;

[n1,n2] = size(f);

r = 11;  % patch size

nt = 4;  % iteration times

% with the above paramter this function takes about 2min.
result = fun_ldmm(f,nt,r); 

subplot(121);imagesc(f);colormap gray;
subplot(122);imagesc(result);colormap gray;



