%%  Reference: 
%Authors: Yongna Jia; Siwei Yu; Lina Liu; Jianwei Ma. 
%A fast rank reduction algorithm for three-dimensional seismic data interpolation. 
% Journal of applied geophysics, 2016, 132: 137-145.
%% Code author: Yongna Jia

%% Acknowledge the following authors and their code (function: EOR1MP)
%Wang Z, Lai M J, Lu Z, et al. Orthogonal rank-one matrix pursuit for low rank matrix completion[J]. 
%SIAM Journal on Scientific Computing, 2015, 37(1): A488-A514.



%% The code is used as a fast algorithm for 3D seismic data interpolation based on Hankel matrix.


%%data: X_3D; sample_matrix:E_3D; 

warning off;
clear;clc;close all;
load X_3D;   load E_3D; rank=3;
[nt,nx ,ny]=size(X_3D);U=find(E_3D);
M_3D=X_3D.*E_3D;
M_f_3D=fft(M_3D);M_f=zeros(nt,nx,ny);
for i=1:nt
    i;
    %%recovery the real data values and the imaginary values
  [Known1,S_r_real]=recover_real_3D(M_f_3D(i,:,:),ny,rank);
  [Known2,S_r_imag]=recover_imag_3D(M_f_3D(i,:,:),ny,rank);
  %% 
  if sum(Known1)==0
      S_r_real=zeros(nx,ny);
  else
  end
   if sum(Known2)==0
      S_r_imag=zeros(nx,ny);
  else
   end
  %%
  S_r=S_r_real+S_r_imag*1i;
  M_f(i,:,:)=S_r;
end
M_r=real(ifft(M_f));
M_r(U)=X_3D(U);
snr=SNR(M_r(:),X_3D(:));


