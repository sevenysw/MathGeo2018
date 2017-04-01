% DDTF3D by YSW
% Interpolation and denoising of high-dimensional seismic data by learning a tight frame
% Author: Siwei Yu
% Reference: Yu S, Ma J, Zhang X, et al. Interpolation and denoising of high-dimensional seismic data by learning a tight frame[J]. Geophysics, 2015, 80(5):V119-V132.
% Email:  jma@hit.edu.cn
% March 22, 2017

clear;clc;close all;

%% 1,******************* load data **************************
% model 1: X_3D_syn, 128x128x128, model 2: elf3D, 256x256x16
% u_: data without noise, u0: adding noise

load 3Dsbb.mat
% load elf3D; 
n = 64;
[n1,n2,n3] = size(D);
u_ = D(101:164,round(n2/2-n/2+1:n2/2+n/2),round(n3/2-n/2+1:n3/2+n/2));
[n1, n2, n3]=size(u_);

u_ = u_/max(abs(u_(:)));

u0 = u_ + 0.2 * randn(n1,n2,n3);

lam = 1.5; thr = 0.8;

%denoise
out=ddtfdenoise3d(u0,lam,thr);

%display
seishow3D(u0);title('noisy data');
seishow3D(out);title('denoised data');
