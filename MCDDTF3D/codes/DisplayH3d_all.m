function DisplayH3d_all(d0)

%Display 3D dictionary

n1 = size(d0,1);
n2 = size(d0,2);
n = round(n1^(1.0/3.0));      
d0(n1,:)=d0(n1,:)*1.00001;
newd = zeros(n+1,n+1,n+1);
% tm =[1 2 5 17 6 18 21 22];  % 4
% tm =[1 2 9 65 10 66 73 74]; % 8
% tm =[1 2 3 4 5 6 7 8]; % 8
% tm =[1 2 4 6 8 9 10 11 24 25 26 28 31 32 41 44] ;%8-di-elf
tm =[1 2 3 4 6 8 9 12 20 21 22 27 32 54 56 58] ;%8-d-X
% tm = [1 2 8 50 9 51 57 58]; % 7
% tm = [1 2 4 5 8 10 11 18 19 23 24 25 26 30 31 32 33 52 53 54]; % 7-2
  tm = 1:64;
for m = 1:1
    figure;
    for i=1:length(tm)

        subplot(8,8,i);
        newd = zeros(n+1,n+1,n+1);
        % DisplayH3d( reshape(d0(:,16/16*(i-1)+1),n,n,n));
%         p = (m-1)*n*n+(i-1)+1;
        p = tm(i) ;
        newd(1:n,1:n,1:n) = reshape(d0(:,p),n,n,n);
        newd(:,:,n+1) = newd(:,:,n);
        newd(:,n+1,:) = newd(:,n,:);
        newd(n+1,:,:) = newd(n,:,:);
        slice(newd,[1,n+1],[1,n+1],[1,n+1]);%colormap gray;
        
        h=findobj(gca,'linestyle','-');
        h(1)=h(2);
        set(h,'linestyle','none');
%         if (i==1) 
%             alpha(1);
%         else
%             alpha(0.7);
%         end

        axis equal;
        axis tight;

        xlabel('x');
        ylabel('y');
        zlabel('z');
        axis off;
        grid off; 
    end
end



end