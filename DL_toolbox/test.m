%% Simultaneously dictionary learning and denoising for seismic data
% Author: F. Bossmann
% Reference:
% S. Beckouche, J. Ma*, Simultaneously dictionary learning and denoising for seismic data, Geophysics, 2014, 79 (3), A27-A31.
% Email:jma@hit.edu.cn
% April 11st, 2017

%%
clf reset
clc
clear all;
close all;
echo off
options.sparse_coding = 'omp_err';
options.manual=1;
%options.dico_sigm = sigma;
options.centerize = 0;
options.centerize_den = 0;
options.q=1;
options.linearis=0
%options.w=8;%defaul 9

%load X_elf3D;
%M0=X_elf3D(1:256,1:256); clear X_elf3D;
load X1
M0=X1(101:356,101:356);
M0=M0/max(M0(:));
%level=0.2;
%noisy3=level*randn(size(M0));
%save noisy3.mat noisy3
load noisy3
M= M0 + noisy3;
%options.sigm2=MAD(M(:));
SNR_noisy=SNR(M,M0)
figure, imagesc(M0), colormap(gray);
figure, imagesc(M), colormap(gray);

%tic
%u_new=nonlocalTV(M);
%tic
%SNR_NTV=SNR(u_new,M0)
%figure, imagesc(u_new), colormap(gray);


%dt=0.0028;time=30;%0.0004
%tic
%atv_M = atv(M,dt,time);
%toc
%figure, imagesc(atv_M), colormap(gray);
%SNR_atv=SNR(atv_M,M0)


% [C,Ct]=curveletTh(M,0.13);%0.035 for level 0.05; 0.055 for 0.08;0.14 for X1 (0.2)
% figure, imagesc(C), colormap(gray);
% SNR_curvelet=SNR(C,M0)


% W=waveletTh2(M,0.13);%matlab toolbox
% figure, imagesc(W), colormap(gray);
% SNR_wavelet=SNR(W,M0)

%[Yw,C1]=waveletTh(X1,cutoff)

for i=1:2
 tic
[D1,X] = perform_dictionary_learning(M0,options);
%load D1;
%figure, imagesc(D1), colormap(gray);

M = perform_dictionary_denoising(M,D1,options);

toc
figure, imagesc(M), colormap(gray);
SNR_denoised=SNR(M,M0)
M0=M;
%%%%%%%%%%%%%%%%%%%%%%%%
[n1,n2]=size(D1)
H0=zeros(3,12);
L0=zeros(9,3);
B=[];
k=0;
for j=1:10:160
    k=k+1;
    A=reshape(D1(:,j),sqrt(n1),sqrt(n1));A=[A,L0];A=[A;H0];
    B=[B,A];
end
 C=[B(:,1:48);B(:,49:96);B(:,97:144);B(:,145:192)];
figure, imagesc(C(1:45,1:45)), colormap(gray);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n1,n2]=size(D1)
H0=zeros(3,12);
L0=zeros(9,3);
B=[];
k=0;
for j=1:6:150
    k=k+1;
    A=reshape(D1(:,j),sqrt(n1),sqrt(n1));A=[A,L0];A=[A;H0];
    B=[B,A];
end
 C=[B(:,1:60);B(:,61:120);B(:,121:180);B(:,181:240);B(:,241:300)];
figure, imagesc(C(1:57,1:57)), colormap(gray);
end