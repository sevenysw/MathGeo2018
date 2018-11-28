% GMD-F denoise test on linear events
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
%

%% preparations

close all;
clc;
clear all;

load line3.mat;
f0 = d;
f = f0 + 0.15*randn(size(f0));
snr(f0,f)
[n1,n2] = size(f0);

% parameters:
alpha = 1500;       % bandwidth constraint
tau = 0.0;         % Lagrangian multipliers dual ascent time step
% tau = 0.01;
K = 3;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 2;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code

tic;[u, u_hat, omega] = GMD (f, alpha, tau, K, DC, init, tol);toc;

snr(f0,sum(u,3))

%% Visualization

dt = 0.004;
x = 1:n2;
y = (1:n1)*dt;

amin = min(f0(:));
amax = max(f0(:));

figure('Name', 'Input image');
imagesc(x,y,f,[amin,amax]);
colormap gray;colorbar;
xlabel('Trace number','FontSize',12);ylabel('Time (s)','FontSize',12);set(gca,'FontSize',12);


figure('Name', 'Reconstructed composite');
imagesc(x,y,sum(u,3),[amin,amax]);
colormap gray;colorbar;
xlabel('Trace number','FontSize',12);ylabel('Time (s)','FontSize',12);set(gca,'FontSize',12);

figure('Name', 'Noise');
imagesc(x,y,sum(u,3)-f,[amin,amax]);
colormap gray;colorbar;
xlabel('Trace number','FontSize',12);ylabel('Time (s)','FontSize',12);set(gca,'FontSize',12);
figure,seishowfreq(sum(u,3),1,1);

if (1)
figure('Name', 'Input spectrum');
imagesc(abs(fftshift(fft2(f))));
colormap gray;
axis equal;
axis off;
hold on;set(gca,'FontSize',12);
colors = 'rbycmgrbycmg';

for k = 1:size(omega,3)
    plot( ((size(f,2)-1)*(0.5+omega(:,1,k))+1.5), ((size(f,1)-1)*(0.5+omega(:,2,k))+1.5), colors(k) );
    plot( ((size(f,2)-1)*(0.5+omega(end,1,k))+1.5), ((size(f,1)-1)*(0.5+omega(end,2,k))+1.5), [colors(k),'o'] );
  
end


end