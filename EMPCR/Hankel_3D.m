function M_H_nx= Hankel_3D( A )
 A=squeeze( A);
%  figure;colormap(gray);imagesc(A);
 [nx,ny]=size(A);
 r=ny/2;
M_H_ny=[];
for i=1:nx
    M_ny=Hankel_2D(A(i,:),r);
    M_H_ny=[  M_H_ny,M_ny];
end
   M_H_nx=zeros(nx/2*ny/2,(nx/2+1)*(ny/2+1));
   for k=1:nx/2
       M_H_nx(((k-1)*ny/2+1):k*ny/2,:)=M_H_ny(:,(ny/2+1)*(k-1)+1:(ny/2+1)*(k+nx/2));
   end
%    figure;colormap(gray);imagesc(M_H_nx);
end

