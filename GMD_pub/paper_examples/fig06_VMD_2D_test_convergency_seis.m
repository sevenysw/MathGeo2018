% GMD-F convergency test on linear events
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
load line3.mat
f = d;

% parameters:
alpha = 2000;       % bandwidth constraint
tau = 0.025;         % Lagrangian multipliers dual ascent time step
% tau = 0.01;
K = 3;              % number of modes
DC = 0;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code


colors = 'rbcymgrbycmg';
no = 10; % 100 in the paper
myomega = cell(no,1);
for i=1:no
    
    i
    tic
    [u, u_hat, omega] = GMD(f, alpha, tau, K, DC, init, tol);
    toc
    myomega{i} = omega;

end

alpha = zeros(size(omega,1),size(omega,3));
for i=1:no
    omega =   myomega{i} ;
    [~, ind] =sort(    omega(end,1,:),'ascend');
    omegas = omega(:,:,ind);

    for k = 1:size(omega,3)
        semilogy(  omegas(:,1,k) , 1:size(omega,1), colors(k) );hold on;
        semilogy(  omegas(end,1,k) , size(omega,1), [colors(k) 'o'] );hold on;
    end
end
xlabel('w_x');ylabel('Iterations');

