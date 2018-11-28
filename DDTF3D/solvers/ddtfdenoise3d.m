function out=ddtfdenoise3d(in,lam,thr)

%denoise with ddtf
%in:    input data
%lam:   training paramter
%thre:  denoise parameter

%out:   denoised result

lvl = 2;wname = 'db1';
[n1,n2,n3] = size(in);
[d0,r] = gen_dic_by_iwt_3d(2^lvl,wname); % initial dictionary with inverse wavelet transform

tic;d = train3d(d0,in,r,lam);

g = pf3d(in,r);

coef = d'*g;
coef = coef.* (abs(coef)>thr);
f    = d*coef;

out = pb3d(f,n1,n2,n3);

end