function seishow3D(D,c,dmin,dmax)

if (~exist('c' ,'var')) 
    c=100;
end

cm = min(D(:));
[n1,n2,n3] = size(D);

di = 0.05; dc = 0.05; dt = 0.004;
% xt = [1, 50, 100, 129, 178, 228];
% xtl= [0, 1, 2, 0, 1 ,2];
% yt = [1, 50, 100, 129, 178, 228];
% ytl= [ 0, 1, 2, 0, 0.2, 0.4];
xt = [1, 50, 100, 150, 200];
xtl= [[0, 50, 100]*di, [50, 100]*dc];
yt = [1, 50, 100, 150, 200, 250];
ytl= [[0, 50, 100]*di, [50, 100, 150]*dt];

I                   = zeros(n1+n3+1,n2+n3+1);
I(1:n3,1:n2)        = squeeze( D(round(n1/2),:,:) )';
I(:,n2+1)  = cm;
I(1:n3,n2+2:n2+n3+1)  = mean(D(:))*ones(n3,n3);
I(n3+2:n3+n1+1,1:n2)  = squeeze(D(:,:, round(n3/2)));
I(n3+1,:)  = cm;
I(n3+2:n3+n1+1,n2+2:n2+n3+1) = squeeze( D(:,round(n2/2),:));
%figure,
I = clip(I,c,c);

if (~exist('dmin' ,'var')) 
    dmin = min(I(:)) ;
end
if (~exist('dmax' ,'var')) 
    dmax = max(I(:)) ;
end

figure,imagesc(I,[dmin,dmax]);axis image;colorbar;  colormap gray;
% set(gca,'XTick',xt,'YTick',yt);
% set(gca,'XTickLabel',xtl,'YTickLabel',ytl);
set(gca,'FontSize',12);
xlabel('Midpoint (km)                   Offset (km)');
ylabel('Time (s)                          Offset (km)');
%imagesc(I);axis image; colormap gray; axis off;colorbar;
% figure, plot(I(:,64));
end
