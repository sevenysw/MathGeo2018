function[P] = extract_patch(M,w,m,pos)

% M : Image to extract patches from
% w : patch size
% m : number of patches to extract
% pos : optional patch position list
% Random position if not specified as input
% Simon Beckouche
% simon@beckouche.fr

[n1,n2] = size(M);

if nargin < 4
    pos = randperm(n1*n2);
    pos = pos(1:m);
end


% Extract patches
P = zeros(w,w,m);
for k = 1:m
    cx = mod(pos(k)-1,n1)+1;
    cy = ceil(pos(k)/n1);
    
    x = cx:cx+w-1;
    y = cy:cy+w-1;
    
    x(x>n1) = 2*n1-x(x>n1);
    y(y>n2) = 2*n2-y(y>n2);
    
    
    P(:,:,k) = M(x,y);
    
    
end