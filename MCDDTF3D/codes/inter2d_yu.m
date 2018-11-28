% interpolation by ddtf
function out= inter2d_yu(out1,raw,mask,d,lambda)
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
[n1,n2]= size(raw);
iter = 20;
r = round((size(d,1))^(1/2));
kk = 1;
%interpolation with iteration threshhold method and low frequency constrain.
for it = 1:iter
it

u_ = raw * 0;
w  = raw * 0;

g = pf2d(out1,r);
coef = d'*g;
coef = coef.* (abs(coef)>(lambda*kk));
f    = d*coef;
u_ = u_ + pb2d(f,n1,n2);
w  = w  + ones(n1,n2);

u_ = u_./w;
kk = kk * 0.9;

out1 = kk*(raw-mask.*u_)+u_;  % with noise
% out1 = (raw-mask.*u_)+u_;       % without noise
end

out = u_;

end