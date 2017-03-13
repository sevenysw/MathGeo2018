function u = pb3d(p,n1,n2,n3) 

%backward patch transform
%p:     patch
%n1~n3: data dimensino
%u:     output data

u_ = zeros(n1,n2,n3);uw = u_;
r = round(size(p,1)^(1/3));
l = 1;
s = 2;
for i=1:s:n1-r+1
    for j=1:s:n2-r+1
        for k=1:s:n3-r+1
            u_(i:i+r-1,j:j+r-1,k:k+r-1)=u_(i:i+r-1,j:j+r-1,k:k+r-1)+reshape(p(:,l),r,r,r);
            uw(i:i+r-1,j:j+r-1,k:k+r-1)=uw(i:i+r-1,j:j+r-1,k:k+r-1)+1;
            l = l + 1;
        end
    end
end

u = u_./uw;

end