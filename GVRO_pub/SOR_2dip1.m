function [ out] = SOR_2dip1( data,  Gx, Gt ,pinv_A1,pinv_A2,lambda1,lambda2,m,n,max_iter)
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
alf=1.2;
d1=0.5*data;d2=0.5*data;
% p1=0; alf = 0.1; increment = 1;
Z1=[Gx*d1 Gt*d1]';[X1,Y1,~] = qr(Z1,0);X1=X1(:,1);Y1=Y1(1,:);
Z2=[Gx*d2 Gt*d2]';[X2,Y2,~] = qr(Z2,0);X2=X2(:,1);Y2=Y2(1,:);
% Z1=[ones(size(Gx*data)) zeros(size(Gx*data))]';[r,c]=size(Z1);
% Z2=Z1;
%  X1 = zeros(r,1);   Y1 = eye(1,c); X2 = X1; Y2 = Y1;
%  A=[lambda*Gx;lambda*Gt;eye(length(data))];
for iter = 1:max_iter
  
    Z1=[Gx*d1 Gt*d1]';
    X01 = X1; Y01 = Y1; d01=d1;
    X1 = Z1*Y1';%
    X1 = X1/sqrt(X1'*X1);
    X1=alf*X1+(1-alf)*X01;

    Y1 = X1'*Z1;   
    Y1=alf*Y1+(1-alf)*Y01;

    d1=pinv_A1*[lambda1*X1(1)*Y1';lambda1*X1(2)*Y1';data];
    d1=alf*d1+(1-alf)*d01;
    
        Z2=[Gx*d2 Gt*d2]';
    X02 = X2; Y02 = Y2; d02=d2;
    X2 = Z2*Y2';%
    X2 = X2/sqrt(X2'*X2);
    X2=alf*X2+(1-alf)*X02;

    Y2 = X2'*Z2;   
    Y2=alf*Y2+(1-alf)*Y02;

    d2=pinv_A1*[lambda1*X2(1)*Y2';lambda1*X2(2)*Y2';data-d1];
    d2=alf*d2+(1-alf)*d02;
end %iter
for iter = 1:20
  
    Z1=[Gx*d1 Gt*d1]';
    X01 = X1; Y01 = Y1; d01=d1;
    X1 = Z1*Y1';%
    X1 = X1/sqrt(X1'*X1);
    X1=alf*X1+(1-alf)*X01;

    Y1 = X1'*Z1;   
    Y1=alf*Y1+(1-alf)*Y01;

    d1=pinv_A2*[lambda2*X1(1)*Y1';lambda2*X1(2)*Y1';data];
    d1=alf*d1+(1-alf)*d01;
    
        Z2=[Gx*d2 Gt*d2]';
    X02 = X2; Y02 = Y2; d02=d2;
    X2 = Z2*Y2';%
    X2 = X2/sqrt(X2'*X2);
    X2=alf*X2+(1-alf)*X02;

    Y2 = X2'*Z2;   
    Y2=alf*Y2+(1-alf)*Y02;

    d2=pinv_A2*[lambda2*X2(1)*Y2';lambda2*X2(2)*Y2';data-d1];
    d2=alf*d2+(1-alf)*d02;
end %iter
%      out11=reshape(d1,m,n);  
%      out110=out11(1:end-1,1:end-1);
%      out110=out110(:);
% %      out110=out110-mean(out110);out120=out120-mean(out120);
%      out110=(d1-mean(out110))/sqrt((out110-mean(out110))'*(out110-mean(out110)));
%    
%      G11x=Gx*out110;G11t=Gt*out110;
%      tv1=G11x'*G11x+G11t'*G11t;
judge1=1;judge2=1;
d11=d1/(d1'*d1);
tv1=(Gx*d11)'*(Gx*d11)+ (Gt*d11)'*(Gt*d11);
  energy11=d1'*d1;
  d22=d2/(d2'*d2);
tv2=(Gx*d22)'*(Gx*d22)+ (Gt*d22)'*(Gt*d22);
  energy22=d2'*d2;
  p1=X1(1)/X1(2);p2=X2(1)/X2(2);
     if p1>p2
        pp=p1;
        p1=p2;
        p2=pp;
        pp=d1;
        d1=d2;
        d2=pp;        
    end
 out=[[d1;p1;tv1;energy11;judge1],[ d2;p2;tv2;energy22;judge2]];
end

