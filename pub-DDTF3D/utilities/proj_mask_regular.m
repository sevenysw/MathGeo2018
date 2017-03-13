function mask = proj_mask_regular(m, n, r)
% u, image
% r, data KNOWN ratio
% type: data missing type
%   'r': random missing rows
%   'c': random missing columns
%   'p': random missing pixels

% For ROW COLUMN missing cases, the max gap between two known
% rows/columns is at most 2(1/r-1) with probability r^2

% code duplicate for cases 'r' & 'c'
mask = zeros(m,n);
% p=1/sqrt(r);
p=1/r;

mask(round(1:p:m),:) = 1;

        
end
