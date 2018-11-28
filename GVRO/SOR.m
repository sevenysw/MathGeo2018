function [ out] = SOR( data,  Gx, Gt ,pinv_A,lambda,max_iter)
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
energy=sum(data.^2);
% energy(energy<0.7*mean(energy))=0;
% d1=randn(size(data));
d1=1/2*data;%d2=d1;
d2=1/2*data;
 alf = 0.5;% increment = 1;
Z=[Gx*data Gt*data]';[r,c]=size(Z);datanrm = max(1,norm(data));  
 X1 = zeros(r,1);   Y1 = eye(1,c);   Res = data;   res = datanrm; 
 X2 = zeros(r,1);   Y2 = eye(1,c);
%  A=[lambda*Gx;lambda*Gt;eye(length(data))];
for iter = 1:max_iter
    Z1=[Gx*d1 Gt*d1]';
            X1 = Z1*Y1';
            [X1,~,~] = qr(X1,0);  Y1 = X1'*Z1;
    d1=pinv_A*[lambda*X1(1)*Y1';lambda*X1(2)*Y1';data-d2];
    
    Z2=[Gx*d2 Gt*d2]';
            X2 = Z2*Y2';
            [X2,~,~] = qr(X2,0);  Y2 = X2'*Z2;
    d2=pinv_A*[lambda*X2(1)*Y2';lambda*X2(2)*Y2';data-d1];
end %iter
    
     out11=reshape(d1,m,n);out12=reshape(d2,m,n);     
%      figure;imagesc([[out11;out12], [out11+out12;reshape(data,m,n )]]);axis image;
%      out111=imfilter(out11,1/9*ones(3,3),'replicate');out222=imfilter(out12,1/9*ones(3,3),'replicate');
%      d11=out111(:);d22=out222(:);
%      
% figure;imagesc([[out11;out12], [out11+out12;reshape(data,m,n )]]);axis image;
% figure;imagesc([[out111;out222], [out111+out222;reshape(data,m,n )]]);axis image;
     energy11=d1'*d1;energy22=d2'*d2;
% %      energy1=G11x'*G11x+G11t'*G11t;energy2=G12x'*G12x+G12t'*G12t;
%      energy11=d11'*d11;energy22=d22'*d22;
     out110=out11(1:end-1,1:end-1);out120=out12(1:end-1,1:end-1);
     out110=out110(:);out120=out120(:);
%      out110=out110-mean(out110);out120=out120-mean(out120);
     out110=(d1-mean(out110))/sqrt((out110-mean(out110))'*(out110-mean(out110)));
     out120=(d2-mean(out120))/sqrt((out120-mean(out120))'*(out120-mean(out120)));
%      ppp=reshape(out110,20,20);
%     sum( sum(ppp(1:end-1,1:end-1).^2))
     G11x=Gx*out110;G11t=Gt*out110;G12x=Gx*out120;G12t=Gt*out120; 
% figure;imagesc([reshape(out110,m,n);reshape(out120,m,n);reshape(data,m,n )]);axis image;
%      G111x=Gx*d11;G111t=Gt*d11;G121x=Gx*d22;G122t=Gt*d22; 
%      energy1=G111x'*G111x+G111t'*G111t;energy2=G121x'*G121x+G122t'*G122t;
     tv1=G11x'*G11x+G11t'*G11t;tv2=G12x'*G12x+G12t'*G12t;
%      if (energy1>energy2*2)&&(energy11>energy22)
%          d2=mean(d2)*ones(size(d2));
%      end
%      if (energy2>energy1*2)&&(energy22>energy11)
%          d1=mean(d1)*ones(size(d1));
%      end
%      if (energy11>energy22*4)
%          d2=0*mean(d2)*ones(size(d2));
%      end
%      if (energy22>energy11*4)
%          d1=0*mean(d1)*ones(size(d1));
%      end
judge1=1;judge2=1;
%      if (energy22<250)
%          d2=0*mean(d2)*ones(size(d2));
%          judge2=0;
%      end
%      if (energy11<250)
%          d1=0*mean(d1)*ones(size(d1));
%          judge1=0;
%      end
% 
% if tv1>2
%     d1=0*mean(d1)*ones(size(d1));judge1=0;
% end
% if tv2>(2)
%     d2=0*mean(d2)*ones(size(d2));judge2=0;
% end

% ratio_tv=max(tv_1,tv_2)/min(tv_1,tv_2);
% if tv1/tv2>1.2
%     d1=0*mean(d1)*ones(size(d1));judge1=0;  p1=0;     
% end
% if tv2/tv1>1.2
%     d2=0*mean(d1)*ones(size(d1));judge2=0; p2=0;          
% end
% if tv1/tv2<1.5
%     d1=0*mean(d1)*ones(size(d1));judge1=1;  p1=0;     
% end
% if tv2/tv1<1.5
%     d2=0*mean(d1)*ones(size(d1));judge2=0; p2=0;          
% end
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

