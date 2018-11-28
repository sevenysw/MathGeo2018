function [ image_resto ,p_resto] = write_patches_tv_sor_center( patches, patch_size,r_sub,c_sub )
% "patches" is an input whose columns respond to patches with leftup point
% index stored in "r_sub,c_sub"
%patches=out;energy=out(end-1,:);figure;plot(energy)
tv=patches(end-1,:);%figure(105);plot(tv);pause
judge=(tv~=0);%judge=ones(size(tv));

p=patches(end,:);
patches=patches(1:end-2,:);
% thresh=1.8;
% patches(:,find(abs(tv)>thresh*mean(abs(tv))))=repmat(mean(patches(:,find(abs(tv)>thresh*mean(abs(tv))))),size(patches,1),1);

image_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
D_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
d=ones(size(patches,1),1);
d=reshape(d,patch_size); 
d(:,end)=[];d(end,:)=[];d(:,1)=[];d(1,:)=[];
p_resto=D_resto;
pp=ones(size(patches,1),1);
pp=reshape(pp,patch_size); 
pp(:,end)=[];pp(end,:)=[];pp(:,1)=[];pp(1,:)=[];
%  for j=2:length(c_sub)-2
%      for i=2:length(r_sub)-2
%          if energy(length(r_sub)*(j-1)+i)<0.5*energy_mean
%              index_max=find(energy==max([energy(length(r_sub)*(j-1)+i) energy(length(r_sub)*(j-1)+i-1) energy(length(r_sub)*(j-1)+i+1) ...
%                  energy(length(r_sub)*(j-2)+i-1) energy(length(r_sub)*(j)+i-1)... 
%                  energy(length(r_sub)*(j-2)+i+1) energy(length(r_sub)*(j)+i+1) ...
%                  energy(length(r_sub)*(j-2)+i) energy(length(r_sub)*(j)+i)]));
%              p(length(r_sub)*(j-1)+i)=p(index_max);
%          end
%      end
%  end
 for j=1:length(c_sub)
     for i=1:length(r_sub)
%          d1=reshape(patches(:,length(r_sub)*(j-1)+i),patch_size);
%          d1(:,end)=[];d1(end,:)=[];
%          image_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)=image_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)+d1;
%          D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)=D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)+d;       

%          d1=reshape(patches(:,length(r_sub)*(j-1)+i),patch_size);
         d1=judge(length(r_sub)*(j-1)+i)*reshape(patches(:,length(r_sub)*(j-1)+i),patch_size);
         d1(:,end)=[];d1(end,:)=[];d1(:,1)=[];d1(1,:)=[];
         image_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)=image_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)+d1;
%          D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)=D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)+d;   
         D_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)=D_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)+judge(length(r_sub)*(j-1)+i)*d;   
         p_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)=p_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)+judge(length(r_sub)*(j-1)+i)*p(length(r_sub)*(j-1)+i)*pp;  
         
     end
 end
image_resto=image_resto./(D_resto+eps);p_resto=p_resto./(D_resto+eps);
%  figure;imagesc(D_resto)
end


