function [u, u_hat, omega] = VMDC(signal, alpha, tau, K, DC, init, tol)
% Complex Variational Mode Decomposition
% Authors: Siwei Yu
% Reference: Yu S, Ma J, Zhang X, et al. Interpolation and denoising of high-dimensional seismic data by learning a tight frame[J]. Geophysics, 2015, 80(5):V119-V132.
% Email:  jma@hit.edu.cn
% May 26th, 2017
%
% Input and Parameters:
% ---------------------
% signal  - the time domain signal (1D) to be decomposed
% alpha   - the balancing parameter of the data-fidelity constraint
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%                    1 = all omegas start uniformly distributed
%                    2 = all omegas initialized randomly
% tol     - tolerance of convergence criterion; typically around 1e-6
%
% Output:
% -------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies
%
% When using this code, please do cite our paper:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing (in press)
% please check here for update reference: 
%          http://dx.doi.org/10.1109/TSP.2013.2288675



%---------- Preparations

% Period and sampling frequency of input signal

f = signal;

% Time Domain 0 to T (of mirrored signal)
T = length(f);
fs = 1/T;
t = (1:T)/T;

% Spectral Domain discretization
freqs = t-0.5-1/T;

% Maximum number of iterations (if not converged yet, then it won't anyway)
N = 100;

% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;

% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

% Initialization of omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
             omega_plus(1,i) = (1/(K)*(i-1)-0.5)/2;
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end
% omega_plus(1,1)=-0.5;
% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(N, length(freqs));

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator



% ----------- Main loop for iterative updates




while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(1:T)*(abs(u_hat_plus(n+1, 1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        % center frequencies
        omega_plus(n+1,k) = (freqs(1:T)*(abs(u_hat_plus(n+1, 1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,1:T,k)).^2);
     
    end
    
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter
    n = n+1;
    
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end


%------ Postprocessing and cleanup


% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(T, K);
u_hat((1):T,:) = squeeze(u_hat_plus(N,(1):T,:));

u = zeros(K,length(t));

for k = 1:K
    u(k,:)=ifft(ifftshift(u_hat(:,k)));
end


% recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end

end