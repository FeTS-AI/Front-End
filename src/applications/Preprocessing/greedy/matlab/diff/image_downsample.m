function imnew = image_downsample(image)
% image_downsample: Downsamples image by factor of 2 using FFT
% usage: imnew = image_upsample(image)
%   Image dimensions must be powers of 2 

[m n] = size(image);

F = fftshift(fft2(image));
CF = F(m/4:end-m/4-1,n/4:end-n/4-1);
imnew = abs(ifft2(ifftshift(CF)));
