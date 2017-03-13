function displayH(H)

%display 2d dictionary

%%% try 1
% hsize = sqrt(size(H,1));
% 
% figure, set(gcf, 'units', 'normalized', 'position', [.1 .06 .8 .8]);
% 
% for i = 1:hsize
%     for j = 0:hsize-1
%         
%         subplot(hsize,hsize, i+j*hsize), imagesc(reshape(H(:,i+j*hsize), hsize,hsize)), colormap(gray), axis('image');
%         axis off;
%         
%     end
% end

%%% try 2
hsize = sqrt(size(H,1));
m = hsize^2 + hsize + 1;

M = -0.1*ones(m,m);

for i = 1:hsize
    for j = 0:hsize-1
        
        idx = i + j*hsize;
        h = H(:, idx);
        if (min(h(:))-max(h(:)))==0
            h(:) = h(:)/max(h(:));
        else
            h = h - min(h(:)); h = h/max(h(:));
        end
        h = reshape(h, hsize,hsize);
        
        jj = i + (i-1)*hsize;
        ii = j+1 + j*hsize;
        
        M(ii+1:ii+hsize, jj+1:jj+hsize) = h;
        
    end
end


set(gcf, 'units', 'normalized', 'position', [.3 .2 .4 .6]);
imagesc(M), colormap(gray), axis('image'), axis off;