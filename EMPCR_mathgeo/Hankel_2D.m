function M = Hankel_2D( df1,r )
N=size(df1);
n=N(2);
C=df1(1:r);
R=df1(r:n);
M=hankel(C,R);
end

