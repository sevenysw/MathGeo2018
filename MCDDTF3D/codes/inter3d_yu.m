% interpolation by ddtf
function out= inter3d_yu(out1,raw,mask,d,lambda)
%out1:  initial data
%raw:   subsampled data
%mask:  sampling mask
%d:     dictionary used
%lambda:threshold parameter

%output
%out:   interpolated data

%All materials are copyrighted to HIT University, and are used for 
%academic research purpose only. Unauthorized redistribution of the
%materials is prohibited.
%Author: Siwei Yu
%Date:   Dec,10,2016
[n1,n2,n3]= size(raw);
iter = 20;
r = round((size(d,1))^(1/3));
kk = 1;
%interpolation with iteration threshhold method and low frequency constrain.
for it = 1:iter
it

u_ = raw * 0;
w  = raw * 0;
for i=1:size(out1,1)-r+1

    g = pf3d(out1(i:i+r-1,:,:),r);
    coef = d'*g;
    coef = coef.* (abs(coef)>(lambda*kk));
    f    = d*coef;
    u_(i:i+r-1,:,:) = u_(i:i+r-1,:,:) + pb3d(f,r,n2,n3);
    w(i:i+r-1,:,:)  = w(i:i+r-1,:,:)  + ones(r,n2,n3);
end

u_ = u_./w;
kk = kk * 0.9;

out1 = (raw-mask.*u_)+u_;  % with noise
% out1 = (raw-mask.*u_)+u_;       % without noise
end

out = out1;

end