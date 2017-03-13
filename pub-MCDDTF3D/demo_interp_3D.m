% DDTF3D by YSW 20161210
% 

clear;clc;close all;

%% 1,******************* load data **************************
% model 1: X_3D_syn, 128x128x128, model 2: elf3D, 256x256x16
% u_: data without noise, u0: adding noise

load 3Dsbb.mat
% load X_3D_syn.mat; 
u_ = D(101:164,101:164,101:164)*0.04;
[n1, n2, n3]=size(u_);

u0 = u_ + 0.0 * randn(n1,n2,n3);
%volume_browser(u0);

%% 2,******************* initialize dictionary ************************** 

lvl = 2;wname = 'db1';
[d0,r] = gen_dic_by_iwt_3d(2^lvl,wname); % initial dictionary with inverse wavelet transform

%% 3,******************* sparse sampling ************************** 

ratio = 0.5;
mask = proj_mask(zeros(n2,n3), ratio, 'p');   %irregular sampling
% mask = proj_mask_regular(n2,n3,ratio);  %regular sampling
mask3 = u0; % make 3D sampling matrix
for i=1:n1
    mask3(i,:,:)= mask;
end
u4 = u0.*mask3;

%% 3,******************* pre-interpolation with zero nearest point method ************************** 

out1 = u4;
for i = 1:n1
    temp = u4(i,:,:);
    temp = reshape(temp(:),n2,n3);
    temp = InpaintingInterp2(temp,mask, 'nearest');
    out1(i,:,:) = temp;
end

%% 
%   4,**************** train dictionary with zero order interpolated initial data ******************** 
%   5,**************** interpolation sampled data with trained filters ***********


lambdat =1.5;
d_ = d0;
u5 = out1;
for i=1:1
    tic;d_ = train3d(d_,u5,r,lambdat);
    u5 = inter3d_yu(u5,u4,mask3,d_,0.8);toc;
end


%%   6,**************** plot the results ***********
seishow3D(u4);title('subsampled data');
seishow3D(u5);title('reconstructed data');
 
