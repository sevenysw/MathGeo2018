%% Registration of multi-component seismic waves
% Author: Wang H.
% Reference:
% Wang H, Cheng Y, Ma J. Curvelet-based registration of multi-component seismic waves[J]. Journal of Applied Geophysics, 2014, 104(5):90-96.
% Email:jma@hit.edu.cn

%N：Length of the data, 数据的长度
%alpha, lambda: parameters of the model,模型参数
%t1：step length, 步长
%dpp, dps：compressional wave, shear wave data 纵横波数据
%c：DCT coefficient, 离散余弦系数

% March 22, 2017

function [tal]=multi_match0325
load('dpp.mat');
load('dps.mat');
N=396;
tal0=1:1:396;%initial tal
c=dct(tal0)';%inatial c
ck=c;
tal=idct(c);
ps=zeros(396,1);
ps1=zeros(396,1);
t=1;
t1=0.1;
altha=1;
lambda=1;
iter=850;
T = lambda/(2*1);

%initial S(tal)
for k=2:395
    if tal(k,1)==fix(tal(k,1))
       ps(k,1)=dps(k,1);
    else
       ps(k,1)=dps(k+1,1)*tal(k,1)+dps(k-1,1)*ceil(tal(k,1))-dps(k-1,1)*tal(k,1)-dps(k+1,1)*floor(tal(k,1));
    end
    
end

%initial S1(tal)
for k=2:395
    if tal(k,1)+t1==fix(tal(k,1)+t1)
       ps1(k,1)=dps(k,1);
    else
       ps1(k,1)=dps(k+1,1)*tal(k,1)+dps(k-1,1)*ceil(tal(k,1))-dps(k-1,1)*tal(k,1)-dps(k+1,1)*floor(tal(k,1));
    end
end

%initial A and b
for k=1:396
    B1(k,1)=(ps1(k,1)-ps(k,1))/t1;
    B2(1,k)=cos(pi*(k-0.5)/396);
end
A=B1*B2;
b=dpp-ps+A*c;

%compute c and A by FIST
 for k=1:iter
     tmpx=ck;
     c=wthresh(ck+A'*(b-A*ck)/altha,'s', T);
     tmpt=t;
     t=(1+sqrt(1+4*t^2))/2;
     ck = c + (tmpt-1)/t*(c-tmpx);
     tal=idct(ck);
     
%update S(tal)
     for k=2:395
         if tal(k,1)==fix(tal(k,1))
            ps(k,1)=dps(k,1);
         else
            ps(k,1)=dps(k+1,1)*tal(k,1)+dps(k-1,1)*ceil(tal(k,1))-dps(k-1,1)*tal(k,1)-dps(k+1,1)*floor(tal(k,1));
         end
     
     end
%update S(tal)
     for k=2:395
         if tal(k,1)+t1==fix(tal(k,1)+t1)
            ps1(k,1)=dps(k,1);
         else
            ps1(k,1)=dps(k+1,1)*tal(k,1)+dps(k-1,1)*ceil(tal(k,1))-dps(k-1,1)*tal(k,1)-dps(k+1,1)*floor(tal(k,1));
         end
     end
%update A and b
     for k=1:396
         B1(k,1)=(ps1(k,1)-ps(k,1))/t1;
         B2(1,k)=cos(pi*(k-0.5)/396);
     end
     A=B1*B2;
     b=dpp-ps+A*ck; 
 end

figure(1),plot(dpp);
figure(2),plot(dps);
figure(3),plot(1:396,tal);
figure(4),plot(50:257,dpp(50:257,1));
figure(5),plot(tal,ps);
        
        
        
        
        

