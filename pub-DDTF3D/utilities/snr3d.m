function snr = snr3d(g, f)
% g: ground truth image
% f: noisy/restored image

g = double(g); % in case of data format is unit8,12,16
f = double(f);

    
    ps = norm(f(:)); % average power of signal
    pn = norm(f(:)-g(:)); % average power of noise
    snr = 20.*log10(ps/pn);
    

    
end
