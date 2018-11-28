function [u, u_hat, omega] = GMD(signal, alpha, tau, K, DC, init, tol)
% GMD-F
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
[Hy,Hx] = size(signal);
[X,Y] = meshgrid((1:Hx)/Hx, (1:Hy)/Hy);


% Spectral Domain discretization
fx = 1/Hx;
fy = 1/Hy;
freqs_1 = X - 0.5 - fx;
freqs_2 = Y - 0.5 - fy;

% N is the maximum number of iterations
N=1000;

% For future generalizations: alpha might be individual for each mode
Alpha = alpha*ones(K,1);

% Construct f and f_hat
f_hat = fftshift(fft2(signal));

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
            omega(1,1,k) = 0.25*cos(pi*(k -1)/maxK);
            omega(1,2,k) = 0.25*sin(pi*(k -1)/maxK);
        end
        
        % Case 1: random on half-plane
    case 1
        for k=1:K
            omega(1,1,k) = rand()-1/2;
            omega(1,2,k) = rand()/2;
        end
        
        % DC component (if expected)
        if DC == 1
            omega(1,1,1) = 0;
            omega(1,2,1) = 0;
        end
        
        % Case 2, using mphi on half-plane
    case 2
        [T1,T2] = size(signal);
        fsignal = fftshift(fft(signal),1);
        ty = round(T1/2+T1/6);
        omega(1,2,:) = (ty-1)/(T1-1)-0.5;       
        omega(1,1,:) = (fun_mpfi(fsignal(ty,:),K)-1)/(T2-1) - 0.5;
    case 3
        omega(1,1,:) = 0;
        omega(1,2,:) = 0;
end


alpha = zeros(N,K);
for k=1:K
    alpha(1,k) = angle(-1i*squeeze(omega(1,1,k))+squeeze(omega(1,2,k)));
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
    
%     compute the halfplane mask for the 2D "analytic signal"
%     HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
    HilbertMask = (sign(freqs_1*0 + freqs_2*(1))+1);
    HilbertMask(ceil((end+1)/2),:) = 0;
    % update first mode accumulator
    sum_uk = u_hat(:,:,end) + sum_uk - u_hat(:,:,k);
    
    % update first mode's spectrum through wiener filter (on half plane)
    u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 * cos(alpha(n,k))+freqs_2 *sin(alpha(n,k))).^2));
%     imagesc(abs((f_hat - sum_uk - mu_hat(:,:)/2)));
    % update first mode's central frequency as spectral center of gravity
    if ~DC
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
    end
    alpha(n+1,k) = angle(-1i*squeeze(omega(n+1,1,k))+squeeze(omega(n+1,2,k)));
    
    % recover full spectrum from analytic signal
    u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
    
    % work on other modes
    for k=2:K
        
        % recompute Hilbert mask
%         HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
        HilbertMask = (sign(freqs_1*0 + freqs_2*(1))+1);
        HilbertMask(ceil((end+1)/2),:) = 0;        
        % update accumulator
        sum_uk = u_hat(:,:,k-1) + sum_uk - u_hat(:,:,k);
        
        % update signal spectrum
        u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 * cos(alpha(n,k))+freqs_2 *sin(alpha(n,k))).^2));
        
        % update signal frequencies
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
        alpha(n+1,k) = angle(-1i*squeeze(omega(n+1,1,k))+squeeze(omega(n+1,2,k)));
        % recover full spectrum from analytic signal
        u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
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
u = zeros(Hy,Hx,K);
for k=1:K
    u_hat(1:round(end/6),:,k) = 0;
    u_hat(round(end/6*5):end,:,k) = 0;
    u(:,:,k) = real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))));
end;

% Should the omega-history be returned, or just the final results?
%omega = omega(n,:,:);
omega = omega(1:n,:,:);

end

