% RMD1D decomposition test on hyperbolic events
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
f = d; 
f = f + 0.0*randn(size(f));
% parameters:
alpha = 0.01;       % bandwidth constraint
tau = 0.005;         % Lagrangian multipliers dual ascent time step
K = 2;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code
Np = 2;
[u, u_hat, omega] = HRMD1D (f,Np,alpha, tau, K, DC, init, tol);

seishow(f,sum(u,3));
%% Visualization
if (1)
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

for k = 1:size(omega,3)
    plot(  omega(:,1,k) ,  n1/2*ones(size(omega,1),1) , colors(k) );
    plot( omega(end,1,k)*ones(n1,1) , [1:n1]', [colors(k)] );
  
end

% 
for k=1:size(u,3)
    figure('Name', ['Mode #' num2str(k)]);
    imagesc(u(:,:,k),[amin,amax]);
    colormap gray;
    axis square;
    axis off;
end


end
