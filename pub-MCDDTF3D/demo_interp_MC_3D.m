% MCDDTF3D by YSW
% Monte Carlo data-driven tight frame for seismic data recovery
% Author: Siwei Yu
% Reference: Yu S, Ma J, Osher S. Monte Carlo data-driven tight frame for seismic data recovery[J]. Geophysics, 2016, 81(4):V327-V340.
% Email:  jma@hit.edu.cn
% March 22, 2017

clear;clc;
addpath(genpath('.'));
%% load data;
data_flag = 3;
if (data_flag==1)
    load X_3D_syn;
    n = 64;
    [n1,n2,n3] = size(X_3D_syn);
    d = X_3D_syn(1:n,n2/2-n/2+1:n2/2+n/2,n3/2-n/2+1:n3/2+n/2);
    sigma = 20;
    thre = 0.1;
elseif (data_flag==2)
    load elf3D.mat;
    d = elf3D;
    sigma = 10;
    thre = 0.1;
    load test05MC.mat noisy_img;
    load 09MC_Noisy.mat mask;
elseif (data_flag==3)
    load 3Dsbb.mat
    d = D(101:164,101:164,101:164);
    sigma = 20;
end
d = d - min(d(:)); 
clear_img = round(d/max(d(:))*255);

%% parameters
patchSize 	= 8; 										% patch size
stepSize  	= 1;                       					% overlap step of data   
trainnum	= 40000;									% the number of samples selected for learning
lambda_1  	= 3.4 * sigma;            					% lambda for learning dictionary
% lambda_1  	= 5.0 * sigma;
% lambda_2  	= 3.1 * sigma;            					% lambda for denoising by learned dictionary
lambda_2  	= 2 * sigma;
% for interpolation data elf3D, if no noise, use 1, if with noise, use 2.
% change interpolation function 
opts.nIter	= 30;										% number loop for constructing data-driven filter bank
%opts.A 		= (1/patchSize^1.5)*ones(patchSize^3,1);        % pre-input filters  (must be orthogonal)

noisy_img 	= round(clear_img + 0*sigma*randn(size(clear_img))); 	% add noise
noisy_img(noisy_img > 255) = 255; 
noisy_img(noisy_img < 0)   = 0; 	
% put the image into range [0,255]


%% build mask,sample and pre interpolation
[n1,n2,n3]  = size(clear_img);
ratio       = 0.5;
mask        = proj_mask(zeros(n2,n3), ratio, 'p');   
mask3       = clear_img * 0;
for i=1:n1
    mask3(i,:,:)= mask;
end
raw = clear_img.*mask3;
% zero order
inter_img   = raw * 0;
for i = 1:n1
    temp = raw(i,:,:);
    temp = reshape(temp(:),n2,n3);
    temp = InpaintingInterp2(temp,mask, 'nearest');
    inter_img(i,:,:) = temp;
end

%% mcddtf interpolation
PSNRinput 	= snr(clear_img, inter_img); 				% PSNR of noisy image

patch_opt = 1;
thre = 0.0025;
nt = 4;
perc = zeros(nt,1); PSNRoutput = zeros(nt,1);


im_out_save = cell(nt,1);
for i = 1:nt
    tic
    patchData = [];
    if (patch_opt ==1 )
        
        Data = im2colstep(inter_img,[patchSize patchSize patchSize],[patchSize patchSize patchSize]);       
        Data_var = var(Data./repmat(max(abs(Data)),[size(Data,1) 1])); 
        Data_var_max = max(Data_var);
        
        for k=1:size(noisy_img,1)-patchSize+1
            
            Data = im2colstep(inter_img(k:k+patchSize-1,:,:),[patchSize patchSize patchSize],[stepSize stepSize stepSize]);       
            Data_var = var(Data./repmat(max(abs(Data)),[size(Data,1) 1])); 
            Data_rnd = Data_var_max * rand(size(Data_var));
            Data_flag = find(Data_rnd<Data_var*thre);
            patchData 	= [patchData ,Data(:, Data_flag)];  %mc sample
            
        end
        
        perc(i) = size(patchData,2)/size(Data,2)/(size(noisy_img,1)-patchSize+1);
           
    elseif (patch_opt ==2 )
        
        for k=1:size(noisy_img,1)-patchSize+1
            Data = im2colstep(inter_img(k:k+patchSize-1,:,:),[patchSize patchSize patchSize],[stepSize stepSize stepSize]);
            rperm 	= randperm(size(Data, 2));
            trainnum = round(size(Data, 2)*thre);
            patchData 	= [patchData, Data(:, rperm(1:trainnum))];  %random sample
        end
        
    elseif (patch_opt ==0  )
        for k=1:size(noisy_img,1)-patchSize+1
            Data = im2colstep(noisy_img(k:k+patchSize-1,:,:),[patchSize patchSize patchSize],[stepSize stepSize stepSize]);
            patchData 	= [patchData, Data]; 
        end
    end
    tp(i)=toc;
    %% Learning filter bank from image patches
    tic
    learnt_dict  = filter_learning_3D(patchData, lambda_1, opts);
    tl(i)=toc;
    tic
    im_out 		 = inter3d_yu(inter_img, raw, mask3, learnt_dict, lambda_2);
    ti(i)=toc;
    inter_img    = im_out;
    im_out_save{i}  = im_out;
    PSNRoutput(i) 	 = snr(clear_img, round(im_out));
end

%% Display
snr(clear_img,im_out)
seishow3D(clear_img); title('Original data')
seishow3D(raw); title('Subsampled data')
seishow3D(im_out); title('Reconstructed data')


