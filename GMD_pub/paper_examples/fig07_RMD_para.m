% RMD2D decomposition test on hyperbolic events
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

load para.mat
% f = imresize(d,[128 128]);
f = d; 
% mask = f * 0;
% mask(:,1:4:end) = 1;
% f = mask.*f;
f = f + 0.0*randn(size(f));
% parameters:
alpha = 0.005;       % bandwidth constraint
tau = 0.01;         % Lagrangian multipliers dual ascent time step
K = 3;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code
Np = 2;
[u, u_hat, omega] = HRMD (f,Np,alpha, tau, K, DC, init, tol);

seishow(f,sum(u,3));
%% Visualization
amin = min(f(:));
amax = max(f(:));
figure('Name', 'Input image');
imagesc(f,[amin,amax]);
colormap gray;
axis equal;
axis off;
[n1,n2] = size(f);
p  = linspace(0,2,n1);
h = linspace(0,1,n2);

pf = pradon(f,h,p,Np);
pmin = min(pf(:)); pmax = max(pf(:));
figure('Name', 'Input spectrum');
imagesc(pf,[pmin pmax]);
colormap gray;
axis equal;
axis off;
hold on;
colors = 'rbycmgrbycmg';
% for k = 1:size(omega,3)
%     plot( size(f,2)*(0.5+omega(:,1,k)), size(f,1)*(0.5+omega(:,2,k)), colors(k) );
% end
for k = 1:size(omega,3)
    plot(  omega(:,1,k) ,  omega(:,2,k) , colors(k) );
    plot( omega(end,1,k) ,  omega(end,2,k), [colors(k),'o'] );
  
end

% 
for k=1:size(u,3)
    figure('Name', ['Mode #' num2str(k)]);
    imagesc(u(:,:,k),[amin,amax]);
%     imagesc(pradon(u(:,:,k),h,p,Np),[pmin pmax]);
%     imagesc(abs(fftshift(fft2(sum(u(:,:,k),3)))));
    colormap gray;
    axis equal;
    axis off;
end

figure('Name', 'Reconstructed composite');
imagesc(sum(u ,3)-f ,[amin,amax]);
colormap gray;
axis equal;
axis off;

