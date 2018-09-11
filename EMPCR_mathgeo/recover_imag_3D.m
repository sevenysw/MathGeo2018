function [Known,S_r_imag] = recover_imag_3D( A,ny,rank)
A=imag(A);
M_n=Hankel_3D(A);
[m,n]=size(M_n);
Known = find(M_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%EOR1MP
[~,~,data]=find(M_n);

[U, Theta, V,~,~] = EOR1MP(m, n, rank, Known, data);
M_r=U*diag(Theta)*V';
M_r(Known)=M_n(Known);
S_r_imag=Antiaverage_B_end(M_r,ny);
end





