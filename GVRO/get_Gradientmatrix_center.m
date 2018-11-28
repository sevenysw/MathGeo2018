function [ Gx, Gy ] = get_Gradientmatrix_center( m,n );
%GET_GRADIENTMATRIX : calculate matrix G_1/2*X=gradient with respect to x/y
%dimension(not include the last row and last column)
% ttt=10;
% data1=get_data('noiz.rsf');
% data1=data1(ttt+[1:10],1:25);figure;imagesc(data1);axis image
% 
% gradient_vec1=[vec(data1(2:end,1:end-1)-data1(1:end-1,1:end-1)) ...
%     vec(data1(1:end-1,2:end)-data1(1:end-1,1:end-1))]';
% figure(1);plot(gradient_vec1(1,:),gradient_vec1(2,:),'r*')
% figure;plot(gradient_vec1(1,:)./gradient_vec1(2,:),'r*')
% data1=rand([m, n])
% data1_vec=vec(data1);
% [m n]=size(data1);

Gx=-diag(ones(m*n,1))+diag(ones(m*n-2*m,1),2*m);
Gx(end-2*m+1:end,:)=[];
index=[1 m];
for i=1:n-3
    index=[index [1 m]+i*m];
end
Gx(index,:)=[];

% plot(Gx*data1_vec-vec(data1(2:end-1,3:end)-data1(2:end-1,1:end-2)));

% Gy=-diag(ones(m*n,1))+diag(ones(m*n-1,1),+1);
% Gy(end-m+1:end,:)=[];
% Gy(m*(1:n-1),:)=[];
Gy=-diag(ones(m*n,1))+diag(ones(m*n-2,1),+2);
Gy(1:m,:)=[];Gy(end-m+1:end,:)=[];
index=[m-1 m];
for i=1:n-3
    index=[index [m-1 m]+i*m];
end
Gy(index,:)=[];
% plot(Gy*data1_vec-vec(data1(3:end,2:end-1)-data1(1:end-2,2:end-1)));
% plot((Gy-Gx)*data1_vec);
%figure;imagesc([Gx]);figure;imagesc([Gy])
end

