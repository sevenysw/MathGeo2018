% the function is used to build a tree and learn the dictionary from the
% knots of the tree
%--------------------------------------------------------------------------
function [Lchild, Rchild,Ind1,Ind2] = Build_Tree_Second(y)
[m1,~,~]=size(y);
y1=y(1,:,:);
s=zeros(m1,1);
for i=1:m1
    y2=y(i,:,:);
    s(i)=norm(y1(:,:)-y2(:,:));
end

[s1,I1]=sort(s);
sumall=sum(s1);
for k=1:m1-1
    suma(k)=sum(s1(1:k)); sumb(k)=sumall-suma(k);
    suma(k)=suma(k)/k; sumb(k)=sumb(k)/(m1-k);
    amin(k)=0;
    for r=1:k
       amin(k)= amin(k)+(s1(r)-suma(k))^2;
    end
    for r=k+1:m1
       amin(k) = amin(k)+(s1(r)-sumb(k))^2;
    end
end
amin=amin;
[~,k]=min(amin);
Ind1=I1(1:k); 
Ind2=I1(k+1:m1);
Lchild=y(Ind1,:,:);
Rchild=y(Ind2,:,:);
end

