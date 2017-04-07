%%  Reference: 
%Authors: Yongna Jia£¬Siwei Yu£¬Lina Liu£¬Jianwei Ma. 
%A fast rank reduction algorithm for three-dimensional seismic data interpolation. 
% Journal of applied geophysics, 2016, 132: 137-145.
%% Code author: Yongna Jia

%% Acknowledge the following authors and their code (function: EOR1MP)
%Wang Z, Lai M J, Lu Z, et al. Orthogonal rank-one matrix pursuit for low rank matrix completion[J]. 
%SIAM Journal on Scientific Computing, 2015, 37(1): A488-A514.



%% The code is used as a fast algorithm for 3D seismic data interpolation based on Hankel matrix.


%%data: original_data_128.mat; sample_matrix:E_3D_cell_originaldata_size; 
%%programe:EMPCR_size_originaldata, figure 5 in the reference.

% Matlab 2012b
warning off;
clear;clc;close all;
load original_data_128;   load E_3D_cell_originaldata_size; rank=3;

ntest = 12;
M_r_cell=cell(1,ntest);Time_cycle=zeros(1,ntest);snr_cycle=zeros(1,ntest);
for i_number=1:ntest
    m=i_number*4
    X_3D=d3_128(1:64,1:m,1:m);
    [nt,nx ,ny]=size(X_3D);
    E_3D=E_3D_cell{1,i_number};U=find(E_3D);
    M_3D=X_3D.*E_3D;
    M_f_3D=fft(M_3D);M_f=zeros(nt,nx,ny);
    t0=clock;
for i=1:nt
    i;
    %%recovery the real data values and the imaginary values
  S_r_real=recover_real_3D(M_f_3D(i,:,:),ny,rank);
  S_r_imag=recover_imag_3D(M_f_3D(i,:,:),ny,rank);
  S_r=S_r_real+S_r_imag*1i;
  M_f(i,:,:)=S_r;
end
M_r=real(ifft(M_f));
M_r(U)=X_3D(U);
t1=clock;
t_end=etime(t1,t0);snr=SNR(M_r(:),X_3D(:));
M_r_cell{1,i_number}=M_r;Time_cycle(i_number)=t_end;snr_cycle(i_number)=snr;
end


