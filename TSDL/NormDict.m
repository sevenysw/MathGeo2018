function D = NormDict (A)
    D=zeros(size(A));
    for i=1:size(A,2)
        if norm(A(:,i))~=0,
            D(:,i)=A(:,i)/norm(A(:,i));
        end
    end
end