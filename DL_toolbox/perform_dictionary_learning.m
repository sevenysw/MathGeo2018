function[D,X]=perform_dictionary_learning(Mn,options)
% Learn dictionary from an image or set of training sample
% options.w : patch size, defaut is 9
% options.p : number of atoms in the dictionary, default is 2 * 9^2 = 162
% options.nout : number of global iteration, default is 100
% optiojs.dico_sigm : assumed noise level, used for omp_err and iterative thresholding, default is estimated with MAD
% options.T : number of active atoms allowed for omp reconstruction, default is 3
% options.C : gain for denoising when using omp_err, default is 1.15
% options.centerize : centerize training sample before learning, default is 0
% options.keep : keep vector used to initialize the dictionary in the training set, default is 0
% options.random_seed : randomly initialize the dictionary with selecting random training sample, default is 1
% options.image_update : fix parameter mu so that at training samples are denoised at each iteration using : Vnew = (Vold + mu * D*x)/(1 + mu) where x is the sparse coding of Vold in D
% options.training_set : training set manually input, skyping patch extraction

verbose = getoptions(options,'verbose',1); % w is the patch size, W = w*w
if isfield(options,'W') && not(isfield(options,'w'))
    W = options.W;
    w = sqrt(W);
elseif isfield(options,'w') && not(isfield(options,'W'))
    w = options.w;
    W = w^2;
else
    W = getoptions(options,'W',81);
    w = getoptions(options,'w',sqrt(W)); 
end
p = getoptions(options,'p',2*W); % Dictionary size
nin = getoptions(options,'nin',25); % iteration number for iterative thresholding
nout = getoptions(options,'nout',100); % iteration number for learning
dico_sigm = getoptions(options,'dico_sigm',mad(Mn(:))); % Assumed noise variance
sparse_coding = getoptions(options,'sparse_coding','omp'); % sparse coding method
T = getoptions(options,'T',3); % Sparsity target zhen using OMP
C = getoptions(options,'C',1.15); % Gain level for OMP and OMPerr
m = getoptions(options,'m',2*max(20*p,40*W)); % Number of training patches to use
centerize = getoptions(options,'centerize',0);
linear = getoptions(options,'linear',0);
keep = getoptions(options,'keep',0);
random_seed = getoptions(options,'random_ssed',1);
check_dico = getoptions(options,'check_dico',0);
upd = getoptions(options,'update','mod');
image_update = getoptions(options,'image_update',0);
if image_update
    update_weight = getoptions(options,'update_weight', dico_sigm^2/30);
end
signe = getoptions(options,'signe',0);

if not(isfield(options,'training_set'))
    n=length(Mn);
    [n1,n2] = size(Mn);
    % Training patches location
    if m == 0 || m>=n1*n2
        m = n1*n2;
        pos = 1:length(Mn(:));
        P = extract_patch(Mn,w,m,pos);
    else
        m = m+p;
        pos = randperm(n1*n2);
        pos = pos(1:m);
        P = extract_patch(Mn,w,m,pos);
    end
    P = P - repmat( mean(mean(P)), [w w] );
    P = reshape(P, [W m]);
    
else
    P = options.training_set;
    m = size(P,2);
end
% Removing empty patches
l=1;
while l<=size(P,2)
    if norm(P(:,l),'fro') == 0
        P = P(:,[1:l-1,l+1:end]);
    else l = l+1;
    end
end



if isfield(options,'D0')
    D =options.D0;
else
    % Removing dictionary patches from training set
    m = size(P,2);
    if random_seed
        sel=randperm(m);
        sel=sel(1:p);
    else
        sel = 1:floor(m/p):p;
    end
    D = P(:,sel);
    if not(keep)
        Pc = zeros(W,m-p);
        l = 1;
        for k = 1:m
            if sum(k == sel)==0
                Pc (:,l)=P(:,k);
                l = l+1;
            end
        end
        P = Pc;
        m = size(P,2);
    end
    
    if centerize
        options.m = m;
        if linear
            if mod(w,2)==0
                options.center_pos = (ceil(W/2))*ones(m,1);
            else
                options.center_pos = (ceil(W/2))*ones(m,1);
            end
        else
            options.center_pos = repmat(ceil(w/2)*ones(1,m),[2,1]);
        end
        P = patch_center(P,options);
        display('Training patches centered')
    end
end
D = D ./ repmat( sqrt(sum(D.^2)), [W, 1] );
if strcmp(sparse_coding,'hard_thresh')
    rho = getoptions(options,'rho',1.4);
    mu = 1.9/norm(D)^2;
    error_target=rho*dico_sigm*w;
    lambda = zeros(m,1) + 1.5*dico_sigm;
end
tic
X=zeros(p,m);
display(['Training patches built (',num2str(m),' patches)',', starting learning...'])
Pinit = P;
% Learning
for iout=1:nout
    time=toc;
    %     display(['Pass ',num2str(iout),'/',num2str(nout),'(',num2str(time),'s)'])
    if mod(10*iout,nout) == 0 && verbose
        display(['Processing ',num2str(floor(100*iout/nout)),'%, ',num2str(floor(time)),'/',num2str(floor(time/iout*nout)),' sec'])
           end
    %% Sparse coding
    if strcmp(sparse_coding,'hard_thresh')
        for iin=1:nin
            X = X + mu*D'*(P-D*X);
            X = max( 1 - repmat(lambda(:)', [p 1]) ./ max(abs(X),1e-10), 0 ) .* X;
            %     X = perform_thresholding( X + mu*D'*(P-D*X), lambda*mu );
            if mod(iin,5)==0
                lambda = lambda * error_target ./ sqrt( sum( (P-D*X).^2 ) )';
                %             lambda_list = [lambda_list lambda(j)];
            end
        end
    elseif strcmp(sparse_coding,'soft_thresh')
        for iin=1:nin
            X = perform_thresholding( X + mu*D'*(P-D*X), lambda*mu, 'soft' );
        end
    elseif strcmp(sparse_coding,'omp')
        options2 = options;
        options2.verbose = 0;
        L = T;
        GAMMA = (OMP(D,P,L,options2));
        X = full(GAMMA);
    elseif strcmp('omp_err',sparse_coding)
        options2 = options;
        options2.verbose = 0;
        epsilon = dico_sigm * C;
        GAMMA = OMPerr(D,P,epsilon,options2);
        X = full(GAMMA);
    elseif strcmp('svt_err',sparse_coding)
        D = D/norm(D);
        opt.stop_crit = 'error_level';
        opt.Th = 'HardTh';
        opt.sparse_min = 'HardTh';
        opt.epsilon = dico_sigm * C;
        alpha = norm(D)^2/1.9;
        lambda = 0.1* max(max(abs(D'*P)));
        lambda_dec = 0.9;
        x0 = X;
        Y = P;
        X =  SVT(x0,Y,D,lambda,lambda_dec,alpha,opt);
    elseif strcmp('svt',sparse_coding)
        D = D/norm(D);
        opt.stop_crit = 'n_it';
        opt.svt_nit = 100;
        opt.sparse_min = 'HardTh';
        opt.Th = 'HardTh';
        alpha = norm(D)^2/1.9;
        x0 = X;
        Y = P;
        X =  SVT(x0,Y,D,lambda,alpha,opt);
    elseif strcmp('admm',sparse_coding)
        %         if iout>1
        %             par.D1init = D1;
        %             par.D2init = D2;
        %             par.B1init = B1;
        %             par.B2init = B2;
        %         end
        clear par;
        phi = pinv(D)/norm(pinv(D));
        options.D = D;
        Y = P;
        Xinit = P;
        A = eye(p);
        AT = eye(w*w);
        par.lambda1 = 0.01;
        par.lambda2 = 0.01;
        par.mu1 = 0.1;
        par.mu2 = 0.1;
        par.nglob = 10;
        par.sparse_min = 'omperr';
        par.nuclear_min = 'SoftTh';
        [Xnew,D1,D2,B1,B2] = nuclear_admm(Xinit,Y,phi,par,options);
        X = D1;
    end
    
    %     %% Nuclear minimization
    %     if nuclear
    %         nr = w;
    %         nc = m;
    %         problem_type = 'NNLS';
    % [II,JJ,bb] = find(P);
    % Jcol = compJcol(JJ);
    % Amap = @(X) sparse_multiply(D,X,II,Jcol);
    % ATmap = @(X) sparse_multiply(pinv(D),X,II,Jcol);
    %
    % %         Amap = @(X)sparse_multiply(D,X);
    % %         ATmap = @(Y)sparse_multiply(pinv(D),Y);
    % % %         ATmap = @(Y)Dinv*Y;
    % %         %         [II,JJ,bb] = find(X);
    %         mutarget= 0.7;
    %         rhotarget = 0;
    %         par = [];
    %         [X,iter,ttime,sd,runhist] =  APGL(nr,nc,problem_type,Amap,ATmap,bb,mutarget,rhotarget,par);
    %         figure;imageplot({reshape(P_c(:,1),[9 9]),reshape(D*Y,[9 9]),reshape(D*X,[9 9])})
    %     end
    %% Dictionary Update
    if max(X(:))>0
        if strcmp(upd,'mod')
            D=P*pinv(X);
        %elseif strcmp(upd,'ksvd')
        %    D = perform_ksvd_update(D,P,options);
        end
        l=1;
        while l<=size(D,2)
            if norm(D(:,l),'fro') == 0
                z = randperm(size(P,2));
                z = z(1);
                %                 D = D(:,[1:l-1,l+1:end]);
                D(:,l) = P(:,z)/norm(P(:,z));
                %                 X = X([1:l-1,l+1:end],:);
            else l = l+1;
            end
        end
        
        p = size(D,2);
        D = D ./ repmat( sqrt(sum(D.^2)), [W, 1] );
        mu = 1.9/norm(D)^2;
    else
        display('X is nul, reducing dico_sigma')
        options.dico_sigm = 0.5*options.dico_sigm;
        iout = 1;
    end
    %% Image update
    if image_update
        P = (update_weight*D*X + Pinit)/(1+update_weight);
    end
end
D = D ./ repmat( sqrt(sum(D.^2)), [W, 1] );
for col = 1:size(D,2)
    D(:,col) = D(:,col) / norm(D(:,col));
end
if signe
    D = repmat(sign(mean(D)),[W 1]).*D;
end
display('Learning complete');
