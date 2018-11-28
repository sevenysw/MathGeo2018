% GMD-F decomposition test on linear events
% Authors: Siwei Yu , Jianwei Ma and Stanley Osher
% siweiyu@hit.edu.cn
% http://homepage.hit.edu.cn/pages/siweiyu
% Initial release 2018-10-01 (c) 2018
%
% When using this code, please do cite our papers:
% -----------------------------------------------
% Yu, Siwei, et al. ¡°Geometric Mode Decomposition.¡± Inverse Problems 
% and Imaging, vol. 12, no. 4, 2018, pp. 831¨C852.
%

%% preparations

close all;
clc;
clear all;

load line3.mat;
f = d;
alpha = 2000;       % bandwidth constraint
tau = 0.025;         % Lagrangian multipliers dual ascent time step
K = 3;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code

tic;[u, u_hat, omega] = GMD (f, alpha, tau, K, DC, init, tol);toc;


%% Visualization
amax = max(f(:));
amin = min(f(:));

figure('Name', 'Input image');
imagesc(f,[amin amax]);
colormap gray;
axis equal;
axis off;

figure('Name', 'Input spectrum');
imagesc(abs(fftshift(fft2(f))));
colormap gray;
axis equal;
axis off;
hold on;
colors = 'rbcymgrbycmg';
% for k = 1:size(omega,3)
%     plot( size(f,2)*(0.5+omega(:,1,k)), size(f,1)*(0.5+omega(:,2,k)), colors(k) );
% end
for k = 1:size(omega,3)
    plot( ((size(f,2)-1)*(0.5+omega(:,1,k))+1.5), ((size(f,1)-1)*(0.5+omega(:,2,k))+1.5), colors(k) ,'LineWidth',1 );
    plot( ((size(f,2)-1)*(0.5+omega(end,1,k))+1.5), ((size(f,1)-1)*(0.5+omega(end,2,k))+1.5), [colors(k),'o'],'LineWidth',1 );
  
end



for k=1:size(u,3)
    figure('Name', ['Mode #' num2str(k)]);
%     imagesc(u(:,:,k),[amin amax]);
    imagesc(abs(fftshift(fft2(sum(u(:,:,k),3)))));
    colormap gray;
    axis equal;
    axis off;
end

figure('Name', 'Reconstructed composite');
imagesc(sum(u(:,:,1:end),3));
colormap gray;
axis equal;
axis off;

