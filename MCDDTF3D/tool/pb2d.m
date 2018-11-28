function u = pb2d(p,n1,n2) 

%backward patch transform
%p:     patch
%n1~n2: data dimensino
%u:     output data

u_ = zeros(n1,n2);uw = u_;
r = round(sqrt(size(p,1)));
k = 1;
for i=1:n1-r+1
    for j=1:n2-r+1        
        u_(i:i+r-1,j:j+r-1)=u_(i:i+r-1,j:j+r-1)+reshape(p(:,k),r,r);
        uw(i:i+r-1,j:j+r-1)=uw(i:i+r-1,j:j+r-1)+1;
        k = k + 1;
    end
end

u = u_./uw;

end