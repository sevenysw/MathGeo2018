function[M2,X] = perform_dictionary_denoising(Mn,D,options)
%% Setting parameters
w = getoptions(options,'w',sqrt(size(D,1)));
p = getoptions(options,'p',sqrt(size(D,2)));
% q = getoptions(options,'q',1);
rho = getoptions(options,'rho',1.4);
manual = getoptions(options,'manual',0);
verbose = getoptions(options,'verbose',1);
sparse_coding = getoptions(options,'sparse_coding','hard_thresh');
C = getoptions(options,'C',1.15);
lam_step = getoptions(options,'lam_step',1);
T = getoptions(options,'T',3);
centerize_den = getoptions(options,'centerize_den',1);
linear = getoptions(options,'linear',0);
if manual
    q = input(['Enter overlap dico_parameter [',num2str(floor(w/2)),'] q = '],'s');
    if isempty(q)
        q = floor(w/2);
    else
        q = str2double(q);
    end
    n_dic = input('Enter number of iterations for denoising [100] n_dic = ','s');
    if isempty(n_dic)
        n_dic = 100;
    else
        n_dic = str2double(n_dic);
    end
    sigm2 = input(['Assumed noise variance [mad = ',num2str(mad(Mn(:))),'] sigm =  '], 's');
    if isempty(sigm2)
        sigm2 = mad(Mn(:));
    else
        sigm2 = str2double(sigm2);
    end
else
    q = options.q;
    if strcmp(sparse_coding,'hard_thresh')
        n_dic = options.n_dic;
        sigm2 = options.sigm2;
    elseif strcmp(sparse_coding,'omp_err')
        sigm2 = options.sigm2;
    end
end

%% Building patches
n = length(Mn);


[dY,dX] = meshgrid(0:w-1,0:w-1);
[y,x] = meshgrid(1:q:n-w/2, 1:q:n-w/2);
m = size(x(:),1);
Xp = repmat(dX,[1 1 m]) + repmat( reshape(x(:),[1 1 m]), [w w 1]);
Yp = repmat(dY,[1 1 m]) + repmat( reshape(y(:),[1 1 m]), [w w 1]);
% Ensure boundary conditions (reflexion)
Xp(Xp>n) = 2*n-Xp(Xp>n);
Yp(Yp>n) = 2*n-Yp(Yp>n);

% Extract a large sub-set of regularly sampled patches
P0 = Mn(Xp+(Yp-1)*n);
P0 = reshape(P0, [w^2, m]);
% Removing mean
a = mean(P0);
P0 = P0 - repmat( a, [w^2 1] );

% Centering noisy patches before denoising
if centerize_den
    options.m = m;
    if linear
        options.center_pos = (ceil(w^2/2)-ceil(w/2))*ones(m,1);
    else
        options.center_pos = repmat(ceil(w/2)*ones(1,m),[2,1]);
        options.orig = [ceil(w/2);ceil(w/2)];
    end
    [P_c,centers] = patch_center(P0,options);
        options.center_pos = centers;
    display('Noisy patches centered')
else
    P_c = P0;
end

%% Sparse coding
if strcmp('hard_thresh',sparse_coding)
    % Calibrating error target
    error_target = rho*w*sigm2;
    % Initializing sparse coding
    X = zeros(p,m);
    lambda = zeros(m,1) + 1.5*sigm2;
    % X = max( 1 - repmat(lambda(:)', [p 1]) ./ max(abs(X),1e-10), 0 ) .* X;
    % Iterating
    mu = 1.9/norm(D)^2;% 1.9/norm(D)^2
    tic
    for k=1:n_dic
        time = toc;
        if mod(10*k,n_dic) == 0 && verbose
            display(['Processing ',num2str(floor(100*k/n_dic)),'%, ',num2str(floor(time)),'/',num2str(floor(time/k*n_dic)),' sec'])
        end
        X = X + mu*D'*(P_c-D*X);
        X = max( 1 - repmat(lambda(:)', [p 1]) ./ max(abs(X),1e-10), 0 ) .* X;
        if mod(k,lam_step)==0
            lambda = lambda * error_target ./ sqrt( sum( (P_c-D*X).^2 ) )';
        end
    end
elseif strcmp('omp_err',sparse_coding)
    epsilon = sigm2 * C;
    %         params.Edata = sqrt(prod(blocksize)) * sigma * gain;   % target error for omp
    GAMMA = OMPerr(D,P_c,epsilon,options);
    X = full(GAMMA);
elseif strcmp('omp',sparse_coding)
    options.verbose = 1;
    X = OMP(D,P_c,T);
end

%% Averaging
% Approximated patches
P_d = D*X;
PA = reshape(P_d, [w w m]);
if centerize_den
    if linear
        PA = patch_center(PA,options);
    else
        PA = patch_decenter(PA,options);
    end
    display('Denoised patches decentered')
end

% Insert back the mean
PA = PA - repmat( mean(mean(PA)), [w w] );
PA = PA + reshape(repmat( a, [w^2 1] ), [w w m]);
% To obtain the denoising, we average the value of the approximated patches PA that overlap
W = zeros(n,n);
M2 = zeros(n,n);
for i=1:m
    x = Xp(:,:,i); y = Yp(:,:,i);
    M2(x+(y-1)*n) = M2(x+(y-1)*n) + PA(:,:,i);
    W(x+(y-1)*n) = W(x+(y-1)*n) + 1;
end
M2 = M2 ./ W;