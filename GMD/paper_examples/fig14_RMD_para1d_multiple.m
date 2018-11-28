% RMD1D multiple attenuation test on hyperbolic events
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


[D,H] = readsegy('gom_cdp_nmo.su');

h1 = [H.offset];
dt = H(1).dt/1000/1000;
% We will use a subset of the data

f = D(800:1200,:);
[n1,n2] = size(f);
tau1 = [1:n1]*dt;
% parameters:
alpha = 0.00001;       % bandwidth constraint
tau = 0.01;         % Lagrangian multipliers dual ascent time step
% tau = 0.00;
K = 2;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code
Np = 2;
[u, u_hat, omega] = HRMD1D (f,Np,alpha, tau, K, DC, init, tol);

%% Visualization
p  = linspace(-1,2,n1);
h = linspace(0,1,n2);
pf = pradon(f,h,p,Np);
colors = 'rbycmgrbycmg';

figure,
pimage2(h1,tau1,clip(f,50,50));colormap(seismic(2));
xlabel('Offset [ft]');
ylabel('Time [s]');
set(gcf,'Position',[100,100,500,800])

figure,pimage2(p,tau1,clip(pf,50,50));colormap(seismic(2));
xlabel('Residual moveout [s]');hold on;
ylabel('Time [s]');
set(gcf,'Position',[100,100,500,800])

for k = 1:size(omega,3)
    plot( p(round(omega(end,1,k)*ones(n1,1))) , [1:n1]'*dt, [colors(k)],'LineWidth',2 );
end

figure,pimage2(h1,tau1,clip(u(:,:,1),50,50));colormap(seismic(2));
xlabel('Offset [ft]');
ylabel('Time [s]');
set(gcf,'Position',[100,100,500,800])

figure,pimage2(h1,tau1,clip(u(:,:,2),50,50));colormap(seismic(2));
xlabel('Offset [ft]');
ylabel('Time [s]');
set(gcf,'Position',[100,100,500,800])

figure,pimage2(h1,tau1,clip(f-squeeze(u(:,:,2)),50,50));colormap(seismic(2));
xlabel('Offset [ft]');
ylabel('Time [s]');
set(gcf,'Position',[100,100,500,800])

if (0)
amin = min(f(:));
amax = max(f(:));
figure('Name', 'Input image');
imagesc(f,[amin,amax]);
colormap gray;
axis equal;
axis off;


pmin = min(pf(:)); pmax = max(pf(:));

figure('Name', 'Input spectrum');
imagesc(pf,[pmin pmax]);
colormap gray;
axis equal;
axis off;
hold on;
colors = 'rbycmgrbycmg';

for k = 1:size(omega,3)
    plot( omega(end,1,k)*ones(n1,1) , [1:n1]', [colors(k)] );
  
end

% % 
for k=1:size(u,3)
    figure('Name', ['Mode #' num2str(k)]);
    imagesc(u(:,:,k),[amin,amax]);
    colormap gray;
    axis equal;
    axis off;
end

figure('Name', 'Reconstructed composite');
imagesc(sum(u ,3) ,[amin,amax]);
colormap gray;
axis equal;
axis off;
end
