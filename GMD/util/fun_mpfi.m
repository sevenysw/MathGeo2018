function w= fun_mpfi(dd,K)
    w = zeros(1,K);
    for t=1:K
        
        ddf = fftshift(fft(dd));
        [y,ind]=max(abs(ddf));
        tmp = ddf*0;
        tmp(ind) = ddf(ind);
        ddr = ifft(fftshift(tmp));
        dd = dd-ddr;   
        w(t) = ind;
    end
    
    w = sort(w);
    
end
