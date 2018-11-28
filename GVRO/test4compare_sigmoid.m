clear; close all;clc
% data1=get_data('/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid.rsf');data1=1000*data1;%figure;imagesc(data1)
% save('sigmoid.mat','data1')
load('sigmoid.mat', 'data1')   
rng(4)
noiselevel=4;
% noiz=imfilter(noiselevel*randn(size(data1)),fspecial('gaussian',[2 2],0.1));
% figure;imagesc([noiz noiselevel*randn(size(data1))]);colormap gray;axis image;
% data_noisy=data1+noiz;figure;imagesc(data_noisy);colormap gray;axis image;
data_noisy=data1+noiselevel*randn(size(data1));
figure(1);subplot(121);imagesc(data1);colormap gray;axis image;title('clean')
subplot(122);imagesc(data_noisy);colormap gray;axis image;title('noisy')
% write_data(data_noisy,'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_noisy0.rsf')
endd=size(data_noisy,2)-6;
 snr_initial=10*log10((norm(data1(5:end-5,5:end-5),'fro'))/(norm(data1(5:end-5,5:end-5)-data_noisy(5:end-5,5:end-5),'fro')))
%% fx-decon, 40*40 patches with 80% overlap
dt=0.004;flow=0;fhigh=100; N=1;lf=5;mu=0.01;patch_size1=[40,40];

overlap1=[ 0.8 0.8];
    [ patches, r_sub,c_sub ] = get_patches( data_noisy, patch_size1, overlap1 );
    fun = @(block_struct) fx_decon_patch(block_struct.data,dt,lf,mu,flow,fhigh,patch_size1);
    tic
    patches_resto=blockproc(patches,[size(patches,1),1],fun,'UseParallel',0);
    toc
    [ data_fx_decon_patch ] = write_patches( patches_resto, patch_size1,r_sub,c_sub );
   figure(1);subplot(131);imagesc(data1);colormap gray;axis image;title('clean')
subplot(132);imagesc(data_noisy);colormap gray;axis image;title('noisy')
subplot(133);imagesc(data_fx_decon_patch);colormap gray;axis image;title('fx-decon-patch')
% snr_fx_decon_patch=10*log10((norm(data1,'fro'))/(norm(data1-data_fx_decon_patch,'fro')));
snr_fx_decon_patch1=10*log10((norm(data1(5:end-5,5:endd),'fro'))/(norm(data1(5:end-5,5:endd)-data_fx_decon_patch(5:end-5,5:endd),'fro')))
%% GVRO with 1 dip
data_prefilt=imfilter(data_noisy,fspecial('gaussian',[3 3],1));
% figure(2);imagesc(data_prefilt);axis image;
% snr_prefilt=10*log10((norm(data1(5:end-5,5:end-5),'fro'))/(norm(data1(5:end-5,5:end-5)-data_prefilt(5:end-5,5:end-5),'fro')))
patch_size2=[16,16];overlap2=overlap1;m=patch_size2(1);n=patch_size2(2); 
[ Gx, Gt ] = get_Gradientmatrix_center( m,n);%figure;imagesc([Gx,Gt])
lambda=4;
A=[lambda*Gx;lambda*Gt;eye(m*n)];
pinv_A=pinv(A);
snr_SOR1_center_all=zeros(1,30);
[ patches, r_sub,c_sub ] = get_patches( data_prefilt, patch_size2, overlap2 );
max_iter=30;
out=zeros(size(patches,1)+2,1*size(patches,2));
 h = waitbar(0,'Please wait...');
tic
for i=1:size(patches,2)
    waitbar(i/size(patches,2),h)
    [out(:,i)]=SOR_single_sor(patches(:,i), Gx, Gt ,pinv_A,lambda,max_iter);
end
toc
close(h)  
[ d1_SOR_center ,p_image1_noisy_center] = write_patches_tv_sor_center( out, patch_size2,r_sub,c_sub );

% figure(3);imagesc([data1,image_resto_SOR_center;data_fx_decon_patch,data_noisy]);
% figure(4);imagesc(p_image1_noisy_center,[-2 2]);
% export_fig Fig/snr_vs_iteration_initialized.fig -native -transparent
% export_fig /home/ccc/papers/denoise/mutidip_ccc_SOR/Fig/snr_vs_iteration.pdf -native -transparent
image_resto_SOR_center=d1_SOR_center;
% figure;imagesc([[data1],[image_resto_SOR_center-data_noisy];[data_fx_decon_patch-data_noisy],[image_resto_SOR_center-data1]]);axis image;
%  snr_SOR_center=10*log10((norm(data1,'fro'))/(norm(data1-image_resto_SOR_center,'fro')));
snr_SOR1_center=10*log10((norm(data1(5:end-5,5:endd),'fro'))/(norm(data1(5:end-5,5:endd)-image_resto_SOR_center(5:end-5,5:endd),'fro')))

figure(11);subplot(221);imagesc(data1);colormap gray;axis image;title('clean')
subplot(222);imagesc(data_noisy);colormap gray;axis image;title('noisy')
subplot(223);imagesc(data_fx_decon_patch);colormap gray;axis image;title('fx-decon-patch')
subplot(224);imagesc(image_resto_SOR_center);colormap gray;axis image;title('GVRO')

figure(21);subplot(221);imagesc(data1);colormap gray;axis image;title('clean')
subplot(222);imagesc(data_noisy);colormap gray;axis image;title('noisy')
subplot(223);imagesc(data_fx_decon_patch-data_noisy);colormap gray;axis image;title('removed-fx-decon-patch')
subplot(224);imagesc(image_resto_SOR_center-data_noisy);colormap gray;axis image;title('removed-GVRO')

%% GVRO with 2 dip
data_prefilt=imfilter(data_noisy,fspecial('gaussian',[3 3],1));
figure(2);imagesc(data_prefilt);colormap gray;axis image;
patch_size2=[16,16]
overlap2=overlap1;m=patch_size2(1);n=patch_size2(2); 
[ Gx, Gt ] = get_Gradientmatrix_center( m,n);%figure;imagesc([Gx,Gt])

lambda1=4;
 A1=[lambda1*Gx;lambda1*Gt;eye(m*n)];
pinv_A1=pinv(A1);
lambda2=50;
 A2=[lambda2*Gx;lambda2*Gt;eye(m*n)];
pinv_A2=pinv(A2);
max_iter=30;
    [ patches_prefilt, r_sub,c_sub ] = get_patches( data_prefilt, patch_size2, overlap2 );

 h = waitbar(0,'Please wait...');
    out_prefilt=zeros(size(patches_prefilt,1)+4,2*size(patches_prefilt,2));
    tic
        for i=1:size(patches_prefilt,2)
            waitbar(i/size(patches_prefilt,2),h)
            [out_prefilt(:,(2*i-1):(2*i))]=SOR_2dip1(patches_prefilt(:,i), Gx, Gt ,pinv_A1,pinv_A2,lambda1,lambda2,m,n,max_iter);
        end
    toc
    close(h)  
[ d1_SOR_prefilt ] = write_patches_center(out_prefilt(1:end-4,1:2:end), patch_size2,r_sub,c_sub );
[ d2_SOR_prefilt ] = write_patches_center(out_prefilt(1:end-4,2:2:end), patch_size2,r_sub,c_sub );

figure(1);imagesc([[d1_SOR_prefilt;d2_SOR_prefilt],[d1_SOR_prefilt+d2_SOR_prefilt;data_noisy],[data_fx_decon_patch;data1]]);colormap gray;axis image;
          

image_resto_SOR_prefilt=d1_SOR_prefilt+d2_SOR_prefilt;     
figure(1);subplot(221);imagesc(data1);colormap gray;axis image;title('clean')
subplot(222);imagesc(data_noisy);colormap gray;axis image;title('noisy')
subplot(223);imagesc(data_fx_decon_patch);colormap gray;axis image;title('fx-decon-patch')
subplot(224);imagesc(image_resto_SOR_prefilt);colormap gray;axis image;title('GVRO')

figure(2);subplot(221);imagesc(data1);colormap gray;axis image;title('clean')
subplot(222);imagesc(data_noisy);colormap gray;axis image;title('noisy')
subplot(223);imagesc(data_fx_decon_patch-data_noisy);colormap gray;axis image;title('removed-fx-decon-patch')
subplot(224);imagesc(image_resto_SOR_prefilt-data_noisy);colormap gray;axis image;title('removed-GVRO')
snr_SOR1=10*log10((norm(data1(5:end-5,5:endd),'fro'))/(norm(data1(5:end-5,5:endd)-image_resto_SOR_prefilt(5:end-5,5:endd),'fro')))

%% get the result via apfs from Madagarscar
%  apfs=get_data('/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/apfs_sigmoid.rsf');
% %  snr_initial=10*log10((norm(data1(5:end-5,5:end-5),'fro'))/(norm(data1(5:end-5,5:end-5)-data_noisy(5:end-5,5:end-5),'fro')))
%  snr_apfs=10*log10((norm(data1,'fro'))/(norm(data1-apfs,'fro')));
%   snr_apfs1=10*log10((norm(data1(5:end-5,5:end-5),'fro'))/(norm(data1(5:end-5,5:end-5)-apfs(5:end-5,5:end-5),'fro')))
%  snr_fx_decon_patch1 
%  snr_SOR1 
% figure(31);imagesc([apfs,image_resto_SOR_prefilt,data_fx_decon_patch]);colormap gray;axis image;
%% show the results in Madagarscar
% 
% write_data(data1(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_initial.rsf')
% write_data(data_noisy(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_noisy.rsf')
% write_data(data_prefilt(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_prefilt.rsf')
% write_data(image_resto_SOR_prefilt(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_SOR.rsf')
% write_data(apfs(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_apfs.rsf')
% write_data(data_fx_decon_patch(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/sigmoid_fxdecon.rsf')
% 
% write_data(image_resto_SOR_prefilt(6:end-5,5:endd)-data_noisy(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_SOR.rsf')
% write_data(apfs(6:end-5,5:endd)-data_noisy(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_apfs.rsf')
% write_data(data_fx_decon_patch(6:end-5,5:endd)-data_noisy(6:end-5,5:endd),'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_fxdecon.rsf')
% %%
% removed_SOR=get_data('/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_SOR.rsf');
% removed_apfs=get_data('/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_apfs.rsf');
% removed_fxdecon=get_data('/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_fxdecon.rsf');
% 
% max_denoised=max(max([removed_SOR removed_apfs removed_fxdecon]))
% min_denoised=min(min([removed_SOR removed_apfs removed_fxdecon]))
% removed_SOR(1,1)=max_denoised;removed_SOR(end,end)=min_denoised;
% removed_apfs(1,1)=max_denoised;removed_apfs(end,end)=min_denoised;
% removed_fxdecon(1,1)=max_denoised;removed_fxdecon(end,end)=min_denoised;
% % write_data(zeros(size(image_resto_SOR_3dip(2:end-1,2:end-1)-data_noisy(2:end-1,2:end-1))),'/home/ccc/codetry/denoise/5jlu/txyapf/realccc2_2/removed_SOR_3dip.rsf')
% write_data(removed_SOR,'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_SOR.rsf')
% write_data(removed_apfs,'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_apfs.rsf')
% write_data(removed_fxdecon,'/home/ccc/codetry/denoise/5jlu/txyapf/sigmoidccc1/removed_fxdecon.rsf')
% 
% figure;imagesc([removed_SOR,removed_apfs,removed_fxdecon]);colormap gray;