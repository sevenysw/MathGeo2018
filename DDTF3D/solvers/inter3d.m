% interpolation by ddtf
function out= inter3d(out1,raw,mask,ratio,d,r,lambda)

%out1:  initial data
%raw:   subsampled data
%mask:  sampling mask
%ratio: sampling ratio: no use
%d:     dictionary used
%r:     size of dictioanry
%lambda:threshold parameter

%output
%out:   interpolated data

%All materials are copyrighted to HIT University, and are used for 
%academic research purpose only. Unauthorized redistribution of the
%materials is prohibited.
%Author: Siwei Yu
%Date:   Dec,10,2016

% interpolation by zero order.
[n1,n2,n3]= size(raw);
iter = 20;

kk = 1;
%interpolation with iteration threshhold method and low frequency constrain.
for it = 1:iter

g = pf3d(out1,r);

coef = d'*g;
coef = coef.* (abs(coef)>(lambda*kk));
f    = d*coef;

u_ = pb3d(f,n1,n2,n3);
kk = kk * 0.9; 
out1 = kk*(raw-mask.*u_)+u_;
end

out=out1;

end