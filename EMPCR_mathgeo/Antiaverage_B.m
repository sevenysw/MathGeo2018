function A= Antiaverage_B( Q,I,J ,ny)
i=I;j=J;Num=0;Sum=zeros(ny/2,ny/2+1);
N=size(Q);
N=[N(1)/(ny/2),N(2)/(ny/2+1)];
while i>=1  &&  j<=N(2)
    Sum=Sum+Q((i-1)*(ny/2)+1:i*(ny/2),(j-1)*(ny/2+1)+1:j*(ny/2+1));
    Num=Num+1;
    i=i-1;j=j+1;
end
A=Sum/Num;
end

