function [ Gx, Gt ] = get_Gradientmatrix_forward( m,n );
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
Gx=zeros((m-2)*(n-2),m*n);
for j=2:(m-1)
    for i=2:(n-1)
        gx=zeros(m,n);
        gx(i-1,j)=0.25;gx(i,j)=0.5;gx(i+1,j)=0.25;
        gx(i-1,j+1)=-0.250;gx(i,j+1)=-0.5;gx(i+1,j+1)=-0.250;%figure;imagesc(gx)
        Gx((i-1)+(j-2)*(m-2),:)=gx(:)';
        
    end
end
%figure;imagesc(Gx)

% tt=imfilter(data1,x_mask,'replicate');figure;imagesc(tt)
% plot(Gx*data1_vec-vec(tt(2:end-1,2:end-1)));
Gt=zeros((m-2)*(n-2),m*n);
for j=2:(m-1)
    for i=2:(n-1)
        gt=zeros(m,n);
        gt(i,j-1)=0.250;gt(i,j)=0.5;gt(i,j+1)=0.250;
        gt(i+1,j-1)=-0.250;gt(i+1,j)=-0.5;gt(i+1,j+1)=-0.250;%figure;imagesc(gx)
        Gt((i-1)+(j-2)*(m-2),:)=gt(:)';
        
    end
end
%figure;imagesc(Gt)

% tt=imfilter(data1,y_mask,'replicate');figure;imagesc(tt)
% plot(Gt*data1_vec-vec(tt(2:end-1,2:end-1)));

end

