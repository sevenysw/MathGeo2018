function mask = proj_mask(u, r, type)
% u, image
% r, data KNOWN ratio
% type: data missing type
%   'r': random missing rows
%   'c': random missing columns
%   'p': random missing pixels

% For ROW COLUMN missing cases, the max gap between two known
% rows/columns is at most 2(1/r-1) with probability r^2

% code duplicate for cases 'r' & 'c'

[m,n] = size(u);

mask = zeros(m,n);

switch type
    
    case 'r'
        
        gap = 1/r;
        gap_ = ceil(gap);
        K = floor(m*r);
        
        for i=1:K
            j = floor((i-1)*gap) + randperm(gap_,1);
            mask(j,:) = 1;
        end
        
    case 'c'
        
        gap = 1/r;
        gap_ = ceil(gap);
        K = floor(n*r);
        
        for i=1:K
            j = floor((i-1)*gap) + randperm(gap_,1);
            mask(:,j) = 1;
        end
        
    case 'p'
        
        pix = randperm(m*n);
        r = fix(r*m*n);
        pix = pix(1:r);
        mask(pix) = 1;
        
    otherwise % pixel-wise missing
        
        pix = randperm(m*n);
        r = fix(r*m*n);
        pix = pix(1:r);
        mask(pix) = 1;
        
end
