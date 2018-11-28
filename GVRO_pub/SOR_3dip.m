function [ out] = SOR_3dip( data,  Gx, Gt ,pinv_A,lambda,m,n,max_iter)
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
d1=1/3*data;d2=d1;d3=d1;
% p1=0; alf = 0.1; increment = 1;
Z1=[ones(size(Gx*data)) zeros(size(Gx*data))]';[r,c]=size(Z1);
Z2=Z1;Z3=Z1;
 X1 = zeros(r,1);   Y1 = eye(1,c); X2 = X1; Y2 = Y1; X3 = X1; Y3 = Y1;
%  A=[lambda*Gx;lambda*Gt;eye(length(data))];
for iter = 1:max_iter
    Z1=[Gx*d1 Gt*d1]';
    X1 = Z1*Y1';
    [X1,~,~] = qr(X1,0);  Y1 = X1'*Z1;
    d1=pinv_A*[lambda*X1(1)*Y1';lambda*X1(2)*Y1';data-d2-d3];
    
    Z2=[Gx*d2 Gt*d2]';
    X2 = Z2*Y2';
    [X2,~,~] = qr(X2,0);  Y2 = X2'*Z2;
    d2=pinv_A*[lambda*X2(1)*Y2';lambda*X2(2)*Y2';data-d1-d3];
    
    Z3=[Gx*d3 Gt*d3]';
    X3 = Z3*Y3';
    [X3,~,~] = qr(X3,0);  Y3 = X3'*Z3;
    d3=pinv_A*[lambda*X3(1)*Y3';lambda*X3(2)*Y3';data-d1-d2];
end %iter
    
%      out11=reshape(d1,m,n);  
%      out110=out11(1:end-1,1:end-1);
%      out110=out110(:);
% %      out110=out110-mean(out110);out120=out120-mean(out120);
%      out110=(d1-mean(out110))/sqrt((out110-mean(out110))'*(out110-mean(out110)));
%    
%      G11x=Gx*out110;G11t=Gt*out110;
%      tv1=G11x'*G11x+G11t'*G11t;
judge1=1;judge2=1;judge3=1;
d11=d1/(d1'*d1);
tv1=(Gx*d11)'*(Gx*d11)+ (Gt*d11)'*(Gt*d11);
energy11=d1'*d1;

d22=d2/(d2'*d2);
tv2=(Gx*d22)'*(Gx*d22)+ (Gt*d22)'*(Gt*d22);
energy22=d2'*d2;

d33=d3/(d3'*d3);
tv3=(Gx*d33)'*(Gx*d33)+ (Gt*d33)'*(Gt*d33);
energy33=d3'*d3;
  p1=X1(1)/X1(2);p2=X2(1)/X2(2);p3=X3(1)/X3(2);

p_sort=[p1,p2,p3];
[p_sort,I] = sort(p_sort,'descend');
d_sort=[d1,d2,d3];  
d_sort=d_sort(:,I);
tv_sort=[tv1,tv2,tv3];   
tv_sort=tv_sort(:,I);
energy_sort=[energy11,energy22,energy33];   
energy_sort=energy_sort(:,I);

 out=[[d_sort(:,1);p_sort(1);tv_sort(1);energy_sort(1);judge1],[ d_sort(:,2);p_sort(2);tv_sort(2);energy_sort(2);...
     judge2],[ d_sort(:,3);p_sort(3);tv_sort(3);energy_sort(3);judge3]];
end

