function A= Antiaverage(Q,I,J)
j=J;Num=0;Sum=0;i=I;
N=size(Q);
while i>=1&&j<=N(2)
Sum=Sum+Q(i,j);
Num=Num+1;
i=i-1;j=j+1;
end
A=Sum/Num;
end

