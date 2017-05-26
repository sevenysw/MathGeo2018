% Complex VMD denoising by YSW
% Author: Siwei Yu
% Reference: Yu S, Ma J. Complex Variational Mode Decomposition for
% Slop-preserving Denoising, summited to IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING
% Email:  jma@hit.edu.cn
% March 22, 2017

clear 
load dn;
[n1,n2] = size(dn);

% parameters of CVMD
alpha = 5000;        % moderate bandwidth constraint
tau = 0.0;            % noise-tolerance (no strict fidelity enforcement)
K = 5;              % 3 modes
DC = 0;             % no DC part imposed
init = 0;           % initialize omegas uniformly
tol = 1e-7;

% parameters of window
nl = 100;           %space window length
nt = 128;           %time window length
dd = zeros(n1,n2);  
dw = zeros(n1,n2);
tic
for i=1:round(nl*0.6):(n2-round(0.4*nl)) %0.6 and 0.4 are parameter for window overlapping.
    i
    if (i<=n2-nl+1)
        dnw = dn(:,i:i+nl-1);
        nlt = nl;
    else
        dnw = dn(:,i:n2);
        nlt  = n2 - i+1;
    end
    for j=1:round(nt*0.6):(n1-round(0.4*nt))
    j
    if (j<=n1-nt+1)
        dnwt = dnw(j:j+nt-1,:);
        ntt = nt;
    else
        dnwt = dnw(j:n1,:);
        ntt  = n1 - j+1;
    end
    ff = fft(dnwt);
    u_ = zeros(ntt,nlt);
    for wi=1:ntt/2+1
        fr = ff(wi,:);
        [u, u_hat, omega] = VMDC(fr, alpha, tau, K, DC, init, tol);
        [~, sortIndex] = sort(omega(end,:));
        omega = omega(:,sortIndex);
        u_hat = u_hat(:,sortIndex);
        u = u(sortIndex,:);
        u_(wi,:) = sum(u);
    end
    u_(ceil(ntt/2)+2:ntt,:)=conj(flipud(u_(2:floor(ntt/2),:)));
    dd(j:j+ntt-1,i:i+nlt-1) = dd(j:j+ntt-1,i:i+nlt-1) + real(ifft(u_));
    dw(j:j+ntt-1,i:i+nlt-1) = dw(j:j+ntt-1,i:i+nlt-1) + ones(ntt,nlt);
    end
end
toc
result = dd./dw;
subplot(121);imagesc(dn);colormap gray;
subplot(122);imagesc(result);colormap gray;
