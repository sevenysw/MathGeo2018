function [ out] = SOR_single_sor( data,  Gx, Gt ,pinv_A,lambda,max_iter)
%MDCA : Multi Directional component analysis
% min (mu/2) || data - d1||_2 + lambda||G(d1)-p*q||^2
%

%
% Input:      data      - the observed vectorized noisy image
%             inv_mat   - to speed up the calculation, we store
%                          precalculated inverse matrices
%             m,n       - size of data
%             Gx,Gt     - precalculated differential matrices
%             max_iter     - max iteration numbers 
%     ** default values of opts are given in { }.
%
% Output: d1,d2      - output components data
%         p1,p2             - dips
% data=patches(:,100);
% energy(energy<0.7*mean(energy))=0;
alf=1.5;
% d1=randn(size(data));
d1=data;%d2=d1;
% p1=0; alf = 0.5; increment = 1;
% Z=[ones(size(Gx*data)) zeros(size(Gx*data))]';[r,c]=size(Z);datanrm = max(1,norm(data));  
%  X = zeros(r,1);   Y = eye(1,c);   Res = data;   res = datanrm; 
Z=[Gx*d1 Gt*d1]';[X,Y,~] = qr(Z,0);X=X(:,1);Y=Y(1,:);
for iter = 1:max_iter
    Z=[Gx*d1 Gt*d1]';
    X0 = X; Y0 = Y; d0=d1;
%    
%         Zo = Z; 
            X = Z*Y';%             X=alf*Xo+(1-alf)*X;

%         if est_rank == 1
%             [X,~,~] = qr(X,0);  Y = X'*Z; %Y=alf*Yo+(1-alf)*Y;
            X = X/sqrt(X'*X);
            X=alf*X+(1-alf)*X0;

            Y = X'*Z;   
            Y=alf*Y+(1-alf)*Y0;
   
    d1=pinv_A*[lambda*X(1)*Y';lambda*X(2)*Y';data];
     d1=alf*d1+(1-alf)*d0;

end %iter
    
%      out11=reshape(d1,m,n);  
%      out110=out11(1:end-1,1:end-1);
%      out110=out110(:);
% %      out110=out110-mean(out110);out120=out120-mean(out120);
%      out110=(d1-mean(out110))/sqrt((out110-mean(out110))'*(out110-mean(out110)));
%    
%      G11x=Gx*out110;G11t=Gt*out110;
%      tv1=G11x'*G11x+G11t'*G11t;
d11=d1/(d1'*d1);
tv=(Gx*d11)'*(Gx*d11)+ (Gt*d11)'*(Gt*d11);
  energy=d1'*d1;
 out=[d1;tv/energy;X(1)/X(2)];
%  figure;imagesc([reshape(d1,m,n);reshape(d2,m,n);reshape(data,m,n )]);axis image;
end

