function [ out] = SOR_single_2( data, Gx_forward, Gt_forward, Gx_backward, Gt_backward  ,pinv_A,lambda,m,n,max_iter)
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
% d1=randn(size(data));
d1=data;%d2=d1;
p1=0; alf = 0.1; increment = 1;
Z=[ones(2*size(Gx_forward*data)) zeros(2*size(Gt_forward*data))]';[r,c]=size(Z);datanrm = max(1,norm(data));  
 X = zeros(r,1);   Y = eye(1,c);   Res = data;   res = datanrm; 
%  A=[lambda*Gx;lambda*Gt;eye(length(data))];
for iter = 1:max_iter
%     Z=[Gx*d1 Gt*d1]';
%     Xo = X; Yo = Y; Res0 = Res; res0 = res; alf0x = alf;d0=d1;
%    
%         Zo = Z; 
            X = Z*Y';
%         if est_rank == 1
            [X,~,~] = qr(X,0);  Y = X'*Z;
%         elseif DoQR
%             [X,R  ] = qr(X,0);  Y = X'*Z;
%         else
%             Xt = X'; Y = linsolve(Xt*X,Xt*Z,linopts);
%         end
%     d1 = linsolve(A,[lambda*X(1)*Y';lambda*X(2)*Y';data]);
    d1=pinv_A*[lambda*X(1)*Y';lambda*X(2)*Y';data];
    Z=[Gx_forward*d1 Gt_forward*d1;Gx_backward*d1 Gt_backward*d1]';
%     d1=alf*d0+(1-alf)*d1;
%     if abs(X(1)/X(2))>1
%         X(1)=0;
%     end
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
tv=(Gx_forward*d11)'*(Gx_forward*d11)+ (Gt_forward*d11)'*(Gt_forward*d11);
  energy=d1'*d1;
 out=[d1;tv/energy;X(1)/X(2)];
%  figure;imagesc([reshape(d1,m,n);reshape(d2,m,n);reshape(data,m,n )]);axis image;
end

