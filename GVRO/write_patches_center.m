function [ image_resto ] = write_patches_center( patches, patch_size,r_sub,c_sub )
% "patches" is an input whose columns respond to patches with leftup point
% index stored in "r_sub,c_sub"
image_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
D_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
d=ones(size(patches,1),1);
d=reshape(d,patch_size); 
d(:,end)=[];d(end,:)=[];d(:,1)=[];d(1,:)=[];
 for j=1:length(c_sub)
     for i=1:length(r_sub)
         d1=reshape(patches(:,length(r_sub)*(j-1)+i),patch_size);
         d1(:,end)=[];d1(end,:)=[];d1(:,1)=[];d1(1,:)=[];
%          norm_d1=1;
         norm_d1=(norm(d1)~=0);
%          pause
         image_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)=image_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)+norm_d1*d1;
%          D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)=D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
%              +patch_size(2)-2)+d;   
         D_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)=D_resto(r_sub(i)+1:r_sub(i)+patch_size(1)-2,c_sub(j)+1:c_sub(j)...
             +patch_size(2)-2)+norm_d1*d;   
     end
 end
image_resto=image_resto./(D_resto+eps);
% figure;imagesc(D_resto)
end

