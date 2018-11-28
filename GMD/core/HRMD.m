function [u, u_hat, omega] = HRMD(signal,Np, alpha, tau, K, DC, init, tol)
% GMD-R
% Authors: Siwei Yu , Jianwei Ma and Stanley Osher
% siweiyu@hit.edu.cn
% http://homepage.hit.edu.cn/pages/siweiyu
% Initial release 2018-10-01 (c) 2018
% Input and Parameters:
% ---------------------
% signal     - the space domain signal (2D) to be decomposed
% alpha      - the balancing parameter for data fidelity constraint
% tau        - time-step of dual ascent ( pick 0 for noise-slack )
% K          - the number of modes to be recovered
% DC         - true, if the first mode is put and kept at DC (0-freq)
% init       - 0 = all omegas start at 0
%              1 = all omegas start initialized randomly
% tol        - tolerance of convergence criterion; typically around 1e-7
%
% When using this code, please do cite our papers:
% -----------------------------------------------
% Yu, Siwei, et al. ¡°Geometric Mode Decomposition.¡± Inverse Problems 
% and Imaging, vol. 12, no. 4, 2018, pp. 831¨C852.
%


% Resolution of image


% N is the maximum number of iterations
N=2000;

% For future generalizations: alpha might be individual for each mode
Alpha = alpha*ones(K,1);

[n1,n2] = size(signal);
p  = linspace(0,2,n1);
h = linspace(0,1,n2);
% Construct f and f_hat
f_hat = pradon(signal,h,p,Np) ;


[Hy,Hx] = size(f_hat);
[X,Y] = meshgrid((1:Hx) , (1:Hy) );
freqs_1 = X;
freqs_2 = Y;

% Storage matrices for (Fourier) modes. All iterations are not recorded.
u_hat = zeros(Hy,Hx,K);
u_hat_old = u_hat;
sum_uk = 0;

% Storage matrices for (Fourier) Lagrange multiplier.
mu_hat = zeros(Hy,Hx);

% N iterations at most, 2 spatial coordinates, K clusters
omega = zeros(N, 2, K);

% Initialization of omega_k
switch init
    case 0
        % spread omegas radially
        
        % if DC, keep first mode at 0,0
        if DC
            maxK = K-1;
        else
            maxK = K;
        end
        for k = DC + (1:maxK) 
            omega(1,1,k) = 0.25*cos(pi*(k-1)/maxK);
            omega(1,2,k) = 0.25*sin(pi*(k-1)/maxK);
        end
        
        % Case 1: random on half-plane
    case 1
        for k=1:K
            omega(1,1,k) = Hy*rand();
            omega(1,2,k) = Hx*rand();
        end
        
        % DC component (if expected)
        if DC == 1
            omega(1,1,1) = 0;
            omega(1,2,1) = 0;
        end
end

%% Main loop for iterative updates

% Stopping criteria tolerances
uDiff=tol+eps;
omegaDiff = tol+eps;

% first run
n = 1;

% run until convergence or max number of iterations
while( ( uDiff > tol || omegaDiff > tol ) && n < N )
    
    % first things first
    k = 1;
    
   
    % update first mode accumulator
    sum_uk = u_hat(:,:,end) + sum_uk - u_hat(:,:,k);
    
    % update first mode's spectrum through wiener filter (on half plane)
    u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2) )./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
    
    % update first mode's central frequency as spectral center of gravity
    if ~DC
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
    end
    
    % recover full spectrum from analytic signal
    
    % work on other modes
    for k=2:K
        
        % update accumulator
        sum_uk = u_hat(:,:,k-1) + sum_uk - u_hat(:,:,k);
        
        % update signal spectrum
        u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2))./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
        
        % update signal frequencies
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % recover full spectrum from analytic signal
    end
    
    % Gradient ascent for augmented Lagrangian
    mu_hat(:,:) = mu_hat(:,:) + tau*(sum(u_hat,3) - f_hat);
    
    % increment iteration counter
    n = n+1;
    
    % convergence?
    uDiff = eps;
    omegaDiff = eps;
    
    for k=1:K
        omegaDiff = omegaDiff + sum(sum(abs(omega(n,:,:) - omega(n-1,:,:)).^2));
        uDiff = uDiff + sum(sum(1/(Hx*Hy)*(u_hat(:,:,k)-u_hat_old(:,:,k)).*conj((u_hat(:,:,k)-u_hat_old(:,:,k)))));
    end
    
    uDiff = abs(uDiff);
    
    u_hat_old = u_hat;

end


%% Signal Reconstruction

% Inverse Fourier Transform to compute (spatial) modes
u = zeros(size(signal,1),size(signal,2),K);
for k=1:K
    u(:,:,k) = ipradon(squeeze(u_hat(:,:,k)),h,p,Np);
end;

% Should the omega-history be returned, or just the final results?
%omega = omega(n,:,:);
omega = omega(1:n,:,:);

end

