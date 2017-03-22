% Removal of curtaining effects by a variational model with directional forward differences
% Author: Jan Henrik Fitschen
% Reference: Jan Henrik Fitschen, Jianwei Ma, Sebastian Schuffb, Removal of curtaining effects by a variational model with directional forward differences, Computer Vision and Image Understanding
% Email:  jma@hit.edu.cn
% March 22, 2017


function [ u,s,t ] = curtainDropping( F, lambda, mu, mu2, mu3, iterations )
% Function for removing curtaining effects (Fitschen, Ma, Schuff 2015)
% using PDHGMp algorithm (Chambolle, Pock 2011).
% Parameters:
%   F               image (3D)
%   lambda          reg. param. for differences along stripes (x-dir) of s
%   mu              reg. param. for differences orthogonal to the stripes 
%                   (y- and z-dir) of u
%   mu2             reg. param. for second order differences of u in z-dir
%   mu3             reg. param. for bidirectional TV of laminar part t 
%                   (x- and y-dir)
%   iterations      number of iterations

[nx,ny,nz] = size(F);

% Create difference matrices
Dx = spdiags([-ones(nx,1), ones(nx,1)],[0,1],nx,nx);
Dx(end,end)=0;
Dy = spdiags([-ones(ny,1), ones(ny,1)],[0,1],ny,ny);
Dy(end,end)=0;
Dz = spdiags([-ones(ny,1), ones(ny,1)],[0,1],nz,nz);
Dz(end,end)=0;
DDz = spdiags([-ones(ny,1), 2*ones(ny,1), -ones(ny,1)],[-1,0,1],nz,nz);
DDz(1,1:2)=0;
DDz(end,end-1:end)=0;
Dxt = Dx';
Dyt = Dy';
Dzt = Dz';
DDzt = DDz';

% Initialize variables
b1 = zeros([nx,ny,nz]);
b1old = b1;
b2 = zeros([nx,ny,nz]);
b2old = b2;
b3 = zeros([nx,ny,nz]);
b3old = b3;
b4 = zeros([nx,ny,nz]);
b4old = b4;
b5 = zeros([nx,ny,nz]);
b5old = b5;
b6 = zeros([nx,ny,nz]);
b6old = b6;
u = F;
s = zeros(size(F));
t = zeros(size(F));
tau = 1/5;
sigma = 1/5;
for k=1:iterations
    k
    % Step 1: Update u,s,t
    sa =    cellfun(@(x)  x*Dy,num2cell(b2old,[1 2]),'UniformOutput',false);
    sb =    cellfun(@(x)  x*Dz,cellfun(@(x) squeeze(x), num2cell(b3old,[1 3]),'UniformOutput',false),'UniformOutput',false);
    sc =    cellfun(@(x)  x*DDz,cellfun(@(x) squeeze(x), num2cell(b6old,[1 3]),'UniformOutput',false),'UniformOutput',false);
    if nz==1
        u = (u - tau * sigma * (cat(3,sa{:}) + (cat(2,sb{:})) ));
    else
        u = (u - tau * sigma * (cat(3,sa{:}) + permute(cat(3,sb{:}),[1,3,2]) + permute(cat(3,sc{:}),[1,3,2])));
    end
    temp = u;
    
    sa =    cellfun(@(x) Dxt*x,num2cell(b1old,[1 2]),'UniformOutput',false);
    s = (s - tau * sigma * (cat(3,sa{:})));
    temp = temp + s;
    
    sa =    cellfun(@(x) Dxt*x,num2cell(b4old,[1 2]),'UniformOutput',false);
    sb =    cellfun(@(x)  x*Dy,num2cell(b5old,[1 2]),'UniformOutput',false);    
    t = (t - tau * sigma * (cat(3,sa{:}) + cat(3,sb{:})));
    temp = temp + t - F;    

    u = u - 1/3 * temp;
    s = s - 1/3 * temp;
    t = t - 1/3 * temp;
    
    %projection with constraint u \in [0,1]
    u = u .* (u>=0); % set to 0 if smaller than zero
    u = u .* (u<=1) + double(u>1); % set to one if larger than one
    s = s + 0.5 * (u .* (u<0)) + 0.5 * ((u - 1) .* (u>1));
    t = t + 0.5 * (u .* (u<0)) + 0.5 * ((u - 1) .* (u>1));
    
    % Step 2: Update p1, p2
    sa = cellfun(@(x) Dx*x,num2cell(s,[1 2]),'UniformOutput',false);
    v1 = softShrinkage(b1 + cat(3,sa{:}), lambda/sigma);
    b1old = b1;
    b1 = b1 + cat(3,sa{:}) - v1;
    
    sa = cellfun(@(x)  x*Dyt,num2cell(u,[1 2]),'UniformOutput',false);
    sb = cellfun(@(x)  x*Dzt,cellfun(@(x) squeeze(x), num2cell(u,[1 3]),'UniformOutput',false),'UniformOutput',false);
    v2 = (b2 + cat(3,sa{:}));
    if nz==1
        v3 = (b3); 
    else
        v3 = (b3 + permute(cat(3,sb{:}),[1,3,2]));
    end
    vabs = sqrt(v2.^2+v3.^2);
    vabsshrink = softShrinkage(vabs,mu/sigma);
    v2 = v2 .* (vabsshrink ./ (vabs+1e-9.*(vabs==0)));
    v3 = v3 .* (vabsshrink ./ (vabs+1e-9.*(vabs==0)));
    b2old = b2;
    b3old = b3;
    b2 = b2 + cat(3,sa{:}) - v2;
    if nz==1
        b3 = b3 - v3;
    else
        b3 = b3 + permute(cat(3,sb{:}),[1,3,2]) - v3;
    end
    
    sa = cellfun(@(x) Dx*x,num2cell(t,[1 2]),'UniformOutput',false);
    sb = cellfun(@(x)  x*Dyt,num2cell(t,[1 2]),'UniformOutput',false);
    v4 = (b4 + cat(3,sa{:}));
    v5 = (b5 + cat(3,sb{:}));
    vabs = sqrt(v4.^2+v5.^2);
    vabsshrink = softShrinkage(vabs,mu3/sigma);
    v4 = v4 .* (vabsshrink ./ (vabs+1e-9.*(vabs==0)));
    v5 = v5 .* (vabsshrink ./ (vabs+1e-9.*(vabs==0)));
    b4old = b4;
    b5old = b5;   
    b4 = b4 + cat(3,sa{:}) - v4;    
    b5 = b5 + cat(3,sb{:}) - v5;
    
    if nz==1
        sa = cellfun(@(x)  x*0,cellfun(@(x) squeeze(x), num2cell(u,[1 3]),'UniformOutput',false),'UniformOutput',false);
    else
        sa = cellfun(@(x)  x*DDzt,cellfun(@(x) squeeze(x), num2cell(u,[1 3]),'UniformOutput',false),'UniformOutput',false);
    end
    v6 = softShrinkage(b6 + permute(cat(3,sa{:}),[1,3,2]), mu2/sigma);
    b6old = b6; 
    b6 = b6 + permute(cat(3,sa{:}),[1,3,2]) - v6;     
    
    % Step 3: Extrapolation
    b1old = 2*b1 - b1old;
    b2old = 2*b2 - b2old;
    b3old = 2*b3 - b3old;   
    b4old = 2*b4 - b4old;
    b5old = 2*b5 - b5old; 
    b6old = 2*b6 - b6old;   
end

end

