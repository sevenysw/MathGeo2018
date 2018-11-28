clear all;close all;clc
addpath(genpath('.'));
%% 参数设置
S       	= 20;                            		
P_Size   	= 8; 							             
S_Size  	= 6;                       		          
load seismic.mat;
C_img       = data{2}(1:256,1:256); 
[h, w] 	  	= size(C_img);  						
N_img 	    = C_img + S*randn(h, w); 	
N_img(N_img > 255) = 255; 
N_img(N_img < 0)   = 0; 						 			
PSNR_I=20*log10(255/std2(C_img-N_img));
%% 建立树和学习字典
Data  		= im2colstep(C_img, [P_Size, P_Size], [S_Size, S_Size]);
[m,n]=size(Data);
for i=1:n
    T(i,:,:)=reshape(Data(:,i),P_Size,P_Size);
end
T_D{1}{1}=T;
C = T(1, : ,: );
Mt=size(T,1);
for i=2:Mt
    C= C+ T(i, : ,:);
end
C1(:,:)= C/Mt;
[U,~,V]=svd(C1);
D0=U(:,1)*V(:,1)';
d=reshape(D0,size(D0,1)*size(D0,2),1);
D(:,1)=d;
m1=size(T_D{1},2);
for i=1:m1
    S(i)=size(T_D{1}{i},1);
end
S_C=max(S);
k=1;
kk=1;
S_Con=S_C;
num=9;
tic;
while S_Con>num
    m2=size(T_D{k},2);
    k=k+1;
    h=0;
    for j=1:m2
        kk=kk+1;
        if S_C(j)>num
            [Lchild, Rchild,Ind1,Ind2,C1] = Build_Tree(T_D{k-1}{j}); 
%              [Lchild, Rchild,Ind1,Ind2] = Build_Tree_Second(T_D{k-1}{j}); 
             T_D{k}{h+1}=Lchild;
             S_L=size(T_D{k}{h+1},1);             
             C = T_D{k}{h+1}(1, : ,: );
             for i=2:S_L
                 C= C+ T_D{k}{h+1}(i, : ,:);
             end
             C1(:,:)= C/S_L;
             [U_c1,~,~]=svd(C1*C1');
             [~,S_c1,V_c1]=svd(C1'*C1);
             D0=U_c1(:,1)*S_c1(1,1)*V_c1(:,1)';
%              D0=U_c1(:,1)*V_c1(:,1)';
             T_D{k}{h+2}=Rchild;
             S_R=size(T_D{k}{h+2},1);
             C = T_D{k}{h+2}(1, : ,: );
             for i=2:S_R
                 C= C+ T_D{k}{h+2}(i, : ,:);
             end
             C1(:,:)= C/S_R;
             [U_c1,~,~]=svd(C1*C1');
             [~,S_c1,V_c1]=svd(C1'*C1);
             D1=U_c1(:,1)*S_c1(1,1)*V_c1(:,1)';
%              D1=U_c1(:,1)*V_c1(:,1)';
%              [U,S,V]=svd(C1);
%              D1=S(1,1)*U(:,1)*V(:,1)';
             S_LR(h+1)=S_L;
             S_LR(h+2)=S_R;
             d0=reshape(D0,size(D0,1)*size(D0,2),1);
             d1=reshape(D1,size(D1,1)*size(D1,2),1);
             d=d0-d1;
             D(:,kk)=d;
             I_n{k}{h+1}=Ind1;
             I_n{k}{h+2}=Ind2;
             h=h+2;
        else
            kk=kk-1;
        end
    end
%     for ij=1:h
%     S(ij)=size(T_D{k}{ij},1);
%     end
    S_C=S_LR;
    S_Con=max(S_C(1:h));
end
D=NormDict(D);
%% 去噪
params.x=N_img;
% params.x=S_img;
params.blocksize=P_Size;
params.dict=D;
params.sigma=18;
[O_img,gama] = ompdenoise2(params,5);
toc;
%% 画图
figure;
subplot(1,3,1);imagesc(C_img);colormap(seismic(2));
subplot(1,3,2);imagesc(N_img);colormap(seismic(2));
subplot(1,3,3);imagesc(O_img);colormap(seismic(2));
PSNR_O=20*log10(255/std2(C_img-O_img));
figure;plot(C_img(:,150));
hold on
plot(O_img(:,150),'r--');
error=norm(O_img-C_img)/norm(C_img)