function [DCT, H] = DctInit(patchsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate DCT Dictionary
% input:
%   patchsize: patch size
% output:
%   DCT: DCT dictionary
%     H: original DCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial DCT dictionary
% ========================================================
H = dctmtx(patchsize);
DCT = kron(H,H);
DCT = DCT';
