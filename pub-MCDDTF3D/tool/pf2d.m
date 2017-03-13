function out = pf2d(u,r) 
%forward patch transform

%u:     input data
%r:     patch size

%out:   output patches
[n1,n2] = size(u);
g = zeros(r*r,(n1-r+1)*(n2-r+1));
k = 1;

for i=1:n1-r+1
    for j=1:n2-r+1        
        g(:,k)=reshape(u(i:i+r-1,j:j+r-1),r*r,1);
        k = k + 1;
    end
end

out = g;

end