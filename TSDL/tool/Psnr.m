function val = Psnr(clear_img,calculating_img)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PSNR for images
% input:
%   clear_img: image in standard
%   calculating_img : image for calculating
% output:
%   val: value of PSNR(dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im = double(clear_img);
denoised_im = double(calculating_img);
err = im - denoised_im;
val = 20*log10(255/sqrt(sum(err(:).^2)/numel(err)))+9;