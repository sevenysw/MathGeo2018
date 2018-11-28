function [ patches, r_sub,c_sub ] = get_patches_slides( image, patch_size, r_slide,c_slide)
% "patches" is an output whose columns respond to patches with leftup point
% index stored in "r_sub,c_sub"
% overlap: in percent,default,50%
 [m,n]=size(image);
 r_sub=1: r_slide:m;%the left up point's row subindex of each patch
 r_sub(find(r_sub>=m+1-patch_size(1)))=[];
 r_sub=[r_sub m+1-patch_size(1)];
 c_sub=1: c_slide:n;
 c_sub(find(c_sub>=n+1-patch_size(2)))=[];
 c_sub=[c_sub n+1-patch_size(2)];
 patches=zeros(patch_size(1)*patch_size(2),length(r_sub)*length(c_sub));
nr=length(r_sub);
nc=length(c_sub);
 for j=1:nc
     for i=1:nr
         patches(:,nr*(j-1)+i)=reshape(image(r_sub(i):r_sub(i)+patch_size(1)-1,c_sub(j):c_sub(j)...
             +patch_size(2)-1),[patch_size(1)*patch_size(2),1]);
     end
 end
end

