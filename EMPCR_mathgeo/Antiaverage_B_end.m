function Anti_B= Antiaverage_B_end(Q,ny )
F=size(Q);F=[F(1)/(ny/2),F(2)/(ny/2+1)];
Anti_B=[];
for i=1:F(1)
    A=Antiaverage_B(Q,i,1,ny);
    B=Antiaverage_end(A);
    Anti_B=[Anti_B;B];
end
for j=2:F(2)
     C=Antiaverage_B(Q,F(1),j,ny);
     D=Antiaverage_end(C);
    Anti_B=[Anti_B;D];
end
Anti_B;
end

