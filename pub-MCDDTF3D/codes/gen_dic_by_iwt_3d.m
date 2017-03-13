function [out2,r]=gen_dic_by_iwt_3d(n,wname)
%initial 3d dictionary by idwt

%n: size of dictionary
%wname: the dwt type used
%output:
%out2: initial dictionary
%r:    size of dictionary

%All materials are copyrighted to HIT University, and are used for 
%academic research purpose only. Unauthorized redistribution of the
%materials is prohibited.
%Author: Siwei Yu
%Date:   Dec,10,2016
r = n;
f = zeros(n,n,n);
l = round(log(n)/log(2));
c = wavedec3(f,l,wname); 

pos = [ 1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;];    
n2 = n*n;
out = zeros(n2,n2,n2);
out2 = zeros(n*n*n,n*n*n);
t = 1;
for m=1:l+1
    
    if (m==1)
        c.dec{1} = 0*c.dec{1};
        c.dec{1} = 1;
        x = waverec3(c);
        c.dec{1} = 0*c.dec{1};
        out(1:n,1:n,1:n) = x ;
        out2(:,1) = x(:);
        t=t+1;
    else

        for p=1:7
        q = (m-2)*7+p+1;
        ni = size(c.dec{q},1);
        nj = size(c.dec{q},2);
        nk = size(c.dec{q},3);
        i0 = 2^(m-2)*pos(p,1);j0 = 2^(m-2)*pos(p,2);k0 = 2^(m-2)*pos(p,3);
            for i=1:ni
                for j=1:nj
                    for k=1:nk
                   
                        c.dec{q} = 0*c.dec{q};
                        c.dec{q}(i,j,k)=1;
                        x = waverec3(c);
                        c.dec{q} = 0*c.dec{q};
 
                        ib = (i0+i-1)*n +1 ; ie = ib + n-1 ;
                        jb = (j0+j-1)*n +1 ; je = jb + n-1 ;
                        kb = (k0+k-1)*n +1 ; ke = kb + n-1 ;
                        out(ib:ie,jb:je,kb:ke) = x ;
                        out2(:,t) = x(:);
                        t = t + 1;
                    % convert to block i,j position
                    end
                end
            end
       end
    end

end


end