function out = pf3d(u,r) 
%forward patch transform

%u:     input data
%r:     patch size

%out:   output patches
[n1,n2,n3] = size(u);
g = zeros(r*r*r,(n1-r+1)*(n2-r+1)*(n3-r+1));
l = 1;
s = 1;
for i=1:s:n1-r+1
    for j=1:s:n2-r+1        
        for k=1:s:n3-r+1
            g(:,l)=reshape(u(i:i+r-1,j:j+r-1,k:k+r-1),r*r*r,1);
            l = l + 1;
        end
    end
end

out = g;

end