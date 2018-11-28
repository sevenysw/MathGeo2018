function [ image_resto ] = write_patches( patches, patch_size,r_sub,c_sub )
% "patches" is an input whose columns respond to patches with leftup point
% index stored in "r_sub,c_sub"
image_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
D_resto=zeros(r_sub(end)+patch_size(1)-1,c_sub(end)+patch_size(2)-1);
d=ones(size(patches,1),1);
d=reshape(d,patch_size); 
d(:,end)=[];d(end,:)=[];
 for j=1:length(c_sub)
     for i=1:length(r_sub)
         d1=reshape(patches(:,length(r_sub)*(j-1)+i),patch_size);
         d1(:,end)=[];d1(end,:)=[];
         image_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
             +patch_size(2)-2)=image_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
             +patch_size(2)-2)+d1;
         D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
             +patch_size(2)-2)=D_resto(r_sub(i):r_sub(i)+patch_size(1)-2,c_sub(j):c_sub(j)...
             +patch_size(2)-2)+d;         
     end
 end
image_resto=image_resto./(D_resto+eps);
% figure;imagesc(D_resto)
end

