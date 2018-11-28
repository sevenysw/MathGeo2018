function [d_] = train3d(d0,u,r,lambda)

% ddtf dictionary training (no low frequency constrain)

%d0:    initial dictionary
%u:     training data
%r:     patch size
%lambda:training paramter
%d_:    trained dictionary


%using ddtf algorithm
%1,rerange the data into small patches. patch size r*r*r, patches number
%(n1-r+1)*(n2-r+1)*(n3-r+1),the data size is (r*r*r)*((n1-r+1)*(n2-r+1)*(n3-r+1))
[n1,n2,n3] = size(u);
%pf3d: transform original data to patch data
%pb3d: transform patch data to original data
g = pf3d(u,r);
%2, use the samples from data(data driven) to train the dictionary

a = d0;
iter = 25;
r1 = 0.0;r2 = 0.0; % weight coefficient for low frequency item

for i=1:iter
   
    %fix the dictionary to update the coefficient.
    kk = 1-(i-1)/iter;
%     kk =1;
    v = a' * g;
    v = wthresh(v, 'h', kk* lambda);
    
    %fix the coefficient and update the dictionary.

    gv= g * v';    
    [s, ~, x] = svd(gv);
    a = s * x';

end

d_ = a;

end