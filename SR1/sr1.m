% sr1(A) - shifted rank-1 approximation of data A
%
% This function calculates the shifted rank-1 approximation of a given
% matrix (or tensor) A. The shifted rank-1 approximation is given as
%
% A ~ S_lambda(uv*)
%
% where u and v are vectors, lambda is a shift-vector and S is the
% column-shift operator (shifts the k-th column of a matrix by lambda_k)
%
% B = sr1(A);
% [B,u,v,lambda] = sr1(A);
% [B,u,v,lambda,maxval,iter] = sr1(A);
% [...] = sr1(A,max_iter);
% [...] = sr1(A,max_iter,sim_shifts);
% [...] = sr1(A,max_iter,sim_shifts,dsp);
%
% B = sr1(A) returns the approximation matrix B = S_lambda(uv*). Input A
% must be a numeric matrix or tensor. The shift is always performed along
% the columns of the matrix A. If A is a tensor, its second to last
% dimension will be reshaped into one dimension to calculate the shift
% rank-1 approximation. The returned solution is then reshaped back into
% the input dimensions.
%
% [B,u,v,lambda] = sr1(A) also returns the vectors u, v and lambda.
%
% [B,u,v,lambda,maxval,iter] = sr1(A) also returns informations about the
% approximation quality and number of iterations: maxval is a vector of
% length 5 containing the largest singular vector of the input data (1),
% after starting guess (2), after global approximation (3), after local
% optimization (4) and of |A| (5, upper bound for optimum). iter contains
% the number of global and local iterations performed.
%
% [...] = sr1(A,max_iter) defines the maximum number of iterations
% performed in each step (global and local iterations). Default value is
% 500. Set max_iter to a length 2 vector to use different maximum number of
% iterations in global and local step.

% [...] = sr1(A,max_iter,sim_shifts) sets the number of simulateneously
% shifted columns per global iteration. Default value is 1 (coincides with
% algorithm in the paper). Set a higher value to increase speed but may
% return a worse approximation. Set sim_shifts = inf to shift all columns
% at a time.
%
% [...] = sr1(A,max_iter,sim_shifts,dsp) sets the display status. Default
% is dsp = false. Set disp = true to output more informations during the
% approximation.
%
% Paper available at: https://arxiv.org/abs/1810.01681

function [approx,u,v,lambda,maxval,iter] = sr1(A,mx_iter,sim_shifts,dsp)

%% Checking input

    % check number and type of input argument
    if nargin<1
        error('Not enough input arguments!');
    elseif nargin>4
        error('Too many input arguments!');
    elseif ~isnumeric(A)
        error('Input data must be of type numeric!');
    end        

    % check if mx_iter is a scalar number or length 2 vector.
    if nargin < 2
        mx_iter = 500;
    elseif numel(mx_iter) < 1 || numel(mx_iter)>2
        error('Number of maximum iterations must be a scalar or length 2 vector!');
    elseif ~isnumeric(mx_iter) || any(mx_iter<0)
        error('Number of maximum iterations must be a positive numeric type');
    end
    
    % check if disp is bool
    if nargin < 4
        dsp = false;
    elseif ~islogical(dsp) || numel(dsp)>1
        error('disp must either be "true" or "false"!');
    end
    
    % check if sim_shifts is a scalar number.
    if nargin < 3
        sim_shifts = 1;
    elseif numel(sim_shifts) ~= 1 || ~isnumeric(sim_shifts) || sim_shifts<1
        error('Number of simultaneously shifts must be a positive scalar!');
    elseif sim_shifts == inf
        sim_shifts = numel(A)/size(A,1);
    end

    % return real values if matrix is real (because of Fourier transform
    % small numerical errors will appear in complex part that should be 0)
    if isreal(A)
        rl = true;
    else
        rl = false;
    end
    
    % check if input is tensor and reshape, store original size later on
    if ~ismatrix(A)
        tensor = true;
        t_size = size(A);
        A = A(:,:);
    else
        tensor = false;
    end

    % set tracking parameters
    iter = [0,0];
    maxval = [svds(A,1),0,0,0,0];

%% Start of algorithm
    
    if dsp
        disp('--- Shifted rank-1 approximation ---');
        if tensor
            disp(['Input is tensor -> reshape into ' num2str(size(A)) ' matrix']);
        end
        disp(['Maximum number of iterations: ' num2str(mx_iter)]);
        disp(['Simulteneously shifted columns in global step: ' num2str(sim_shifts)]);
    end

    tic;
    % go to Fourier domain
    A = fft(A);
    
    % get sizes
    [M,N] = size(A);
    
    % starting guess
    [~,lambda] = max(abs(ifft(abs(A).*A)));
    lambda = lambda-1;
    
    % shift phase
    A = A.*exp(((2*pi*1i/M)*(0:(M-1)).')*lambda);

    % global amplitude
    [u_opt,s,~] = svds(abs(A),1);
    maxval(5) = s/sqrt(M);
    maxval(2) = svds(A,1)/sqrt(M);
    
    t = toc;
    if dsp
        disp(['upper bound for largest SV: ' num2str(maxval(5))]);
        disp(['largest SV input data: ' num2str(maxval(1))]);
        disp(['Starting guess calculated after ' num2str(t) ' seconds']);
        disp(['largest SV starting guess: ' num2str(maxval(2))]);
    end
        
    
    % global approximation
    tic;
    u_iter = u_opt;
    while true
        
        % get largest left singular vector
        [u_iter,~,~] = svds(A,1,'largest','LeftStartVector',u_iter);

        % caclulate phase (is amplitude is 0, set phase to 1)
        u_iter = u_iter./abs(u_iter);
        u_iter(isnan(u_iter)) = 1;

        % combine phase and amplitude
        u_iter = u_opt.*u_iter;
        
        % calculate matrix of all possible shift updates
        tmpMtx = abs(ifft(A.*conj(u_iter))).^2;
        
        % get maximal difference to actual shift
        [v,s] = max(tmpMtx-tmpMtx(1,:));
        [~,k] = maxk(v,sim_shifts);
        s = s(k);
        
        % check if converged
        if all(s==1)                % no shift -> local optimum
            break;
        elseif iter(1)>= mx_iter(1) % maximum number of iterations -> stop
            break;
        end
        
        % increase iteration counter
        iter(1) = iter(1)+1;
        
        % set new shift
        lambda(k) = lambda(k)+s-1;
        
        % update matrix
        A(:,k) = A(:,k).*exp(((2*pi*1i/M)*(0:(M-1)).')*(s-1));
    end
    
    % set singular value after global approximation
    maxval(3) = svds(A,1)/sqrt(M);
    
    t = toc;
    if dsp
        disp(['Global approximation termiated after ' num2str(iter(1)) ' iterations (' num2str(t) ' seconds)']);
        disp(['largest SV global approximation: ' num2str(maxval(3))]);
    end
    
    % local optimization
    tic;
    while true
        
        % get largest left singular vector
        [u_iter,~,~] = svds(A,1,'largest','LeftStartVector',u_iter);
        
        % calculate matrix of all possible shift updates
        tmpMtx = abs(ifft(A.*conj(u_iter))).^2;
        
        % get maximal difference to actual shift
%         [~,k] = max(vec(tmpMtx-tmpMtx(1,:)));
        
        tmp = tmpMtx-tmpMtx(1,:);
        [~,k] = max(tmp(:));
        
        [s,k] = ind2sub([M,N],k);
        
        % check if converged
        if s==1                 % no shift -> local optimum
            break;
        elseif iter(2)>=mx_iter(end)  % maximum number of iterations -> stop
            break;
        end
        
        % increase iteration counter
        iter(2) = iter(2)+1;
        
        % set new shift
        lambda(k) = lambda(k)+s-1;
        
        % update matrix
        A(:,k) = A(:,k).*exp((2*pi*1i*(s-1)/M)*(0:(M-1)).');
    end
    
    % calculate best rank-1 approximation
    [u,s,v] = svds(A,1,'largest','LeftStartVector',u_iter);
    maxval(4) = s/sqrt(M);
    
    t = toc;
    if dsp
        disp(['Local optimization termiated after ' num2str(iter(2)) ' iterations (' num2str(t) ' seconds)']);
        disp(['largest SV local optimization: ' num2str(maxval(4))]);
    end
    
    % back to time domain, set to real values if input was real
    u = ifft(u);
    if rl
        u = real(u);
        v = real(v);
        if dsp
            disp('Removed complex part of u and v (input was real)');
        end
    end
    
    % normalize u
    u = u*sqrt(length(u));
    v = v*s/sqrt(length(u));
    
    % perform shift in real domain
    tic;
    approx = u*v';
    for k=1:N
        approx(:,k) = circshift(approx(:,k),lambda(k));
    end
    t = toc;
    if dsp
        disp(['Calculate and shift approximation in time domain (avoids fft errors on values) in ' num2str(t) ' seconds']);
    end
        
    
    % reshape if input is tensor
    if tensor
        approx = reshape(approx,t_size);
        v = reshape(v,t_size(2:end));
        lambda = reshape(lambda,t_size(2:end));
        if dsp
            disp('Reshaped output back into tensor dimensions');
        end
    end
end