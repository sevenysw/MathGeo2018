% the function is used to build a tree and learn the dictionary from the
% knots of the tree
%--------------------------------------------------------------------------
function [Lchild, Rchild,Ind1,Ind2,C1] = Build_Tree(y)
 C = y(1, : ,: );
[m1,m2,m3]=size(y);
% caculate the mean of the training patches--------------------------------
for i=2:m1
    C= C+ y(i, : ,:);
end
C1(:,:)= C/m1;
% caculate the left and right covariance matrix CL,CR----------------------
CL= zeros(m2,m3);
for i=1:m1
    B1(:,:)=y(i,:,:)-C/m1;    
    CL=CL+B1*B1';
end
CL=CL/m1;
CR= zeros(m2,m3);
for i=1:m1
    B1(:,:)=y(i,:,:)-C/m1;    
    CR=CR+B1'*B1;
end
CR=CR/m1;
[~,~,V]=svd(CL);
u1=V(:,1);
[~,~,V]=svd(CR);
v1=V(:,1);
s=zeros(m1,1);
for i=1:m1
    yk(:,:)=y(i,:,:);
    s(i)=u1'*yk*v1;
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

