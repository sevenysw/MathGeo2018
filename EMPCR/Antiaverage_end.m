function de= Antiaverage_end(Z)
F=size(Z);de=[];
for i=1:F(1)
    A=Antiaverage(Z,i,1);
   de=[de A];
end
for j=2:F(2)
    A=Antiaverage(Z,F(1),j);
    de=[de A];
end
end

