%version 2, split the data gradually, not as a whole.

clear;clc;

%% load data and normalization;
data_flag = 1;
load prestack.mat;
cmp3 = cmp3 - min(cmp3(:));
cmp3 = round(cmp3/max(cmp3(:))*255);
[n1,n2] = size(cmp3);
d = zeros(n1,n2);
d(:,1:2:end) = cmp3(:,1:2:end);
d(:,2:2:end) = cmp3(:,1:2:end);
clear_img = d;

%% parameters
patchSize 	= 8; 										% patch size
stepSize  	= 1;                       					% overlap step of data   
trainnum	= 40000;									% the number of samples selected for learning
sigma = 15;
lambda_1  	= 5 * sigma;            					% lambda for learning dictionary
lambda_2  	= 0.5 * sigma; 
opts.nIter	= 30;										% number loop for constructing data-driven filter bank
noisy_img = clear_img;

%% sample and pre-interpolation
mask = zeros(n1,n2);
mask(:,1:2:end) = 1;
raw = clear_img.*mask;
inter_img = InpaintingInterp2(raw,mask, 'nearest');
PSNRinput 	= snr(clear_img, inter_img); 				% PSNR of noisy image

%% mcddtf interpolation
patch_opt = 1;  % patching method. 1: monte carlo, 2: random, 3 regular,4 full
nt = 5;         % time of iteration: traning-interpolation
PSNRoutput = zeros(nt,1);

% perc = 0.0118;
im_out_save = cell(nt,1);
for i = 1:nt
    tic
    patchData = [];
   if (patch_opt ==1 )
        
        %Monte Carlo patching
        Data = im2colstep(inter_img,[patchSize patchSize],[patchSize patchSize]);       
        Data_var = var(Data);
        Data_var_max = max(Data_var);
        
        thre = 0.2;
        Data = im2colstep(inter_img,[patchSize patchSize],[stepSize stepSize]);   
        Data2 = im2colstep(inter_img,[patchSize patchSize],[stepSize stepSize]);
        Data_var = var(Data);
        Data_rnd = Data_var_max * rand(size(Data_var));
        Data_flag = find(Data_rnd<(Data_var*thre));
        patchData 	= Data2(:, Data_flag);    
        perc = size(patchData,2)/size(Data,2);
        
           
    elseif (patch_opt ==2 )
        
        Data = im2colstep(inter_img,[patchSize patchSize],[stepSize stepSize]);
        rperm 	= randperm(size(Data, 2));
        trainnum = round(size(Data, 2)*perc);
        Data_flag = rperm(1:trainnum);
        patchData 	= Data(:,Data_flag );  %random sample


     elseif (patch_opt ==3 )
        

        Data = im2colstep(inter_img,[patchSize patchSize],[stepSize stepSize]);
        Data_flag = round(1:(1/perc):size(Data, 2));
        patchData 	= Data(:, Data_flag);  %regular sample


    elseif (patch_opt ==4  )

        Data = im2colstep(inter_img,[patchSize patchSize],[stepSize stepSize]);
        Data_flag = round(1:size(Data, 2));
        patchData 	= Data(:, Data_flag);  %full

    end
    tp(i)=toc;
    %% Learning filter bank from image patches
    tic
    learnt_dict  = filter_learning_2D(patchData, lambda_1, opts);
    tl(i)=toc;
    tic
    im_out 		 = inter2d_yu(inter_img, raw, mask, learnt_dict, lambda_2);
    ti(i)=toc;
    inter_img    = im_out;
    im_out_save{i}  = im_out;
    PSNRoutput(i) 	 = snr(clear_img, round(im_out));
end

% Display

dx = 0.01; dt = 0.004;
x = (61:80)*dx;
z = (401:560)*dt;

% Subsampled
figure('Position',[100 100 240 600])
cmpt = cmp3(401:560,61:80);
a = cmpt;
a = a - mean(a(:));
trmx= max(abs(a));
amx=mean(trmx);
wigb(cmpt,1,x,z,amx);set(gca,'FontSize',12);
xlabel('Distance (km)');ylabel('Time (s)');

% Reconstruction
figure('Position',[100 100 240 600])
cmpr = round(im_out_save{end}(401:560,61:80));
wigb(cmpr,1,x,z,amx);
set(gca,'FontSize',12);
xlabel('Distance (km)');ylabel('Time (s)');

% Error
x = (61:80)*dx;
z = (401:560)*dt;
figure('Position',[100 100 240 600])
cmpdif = cmpr-cmpt;
wigb(cmpdif,1,x,z,amx);
set(gca,'FontSize',12);
xlabel('Distance (km)');ylabel('Time (s)');


