function imnew = image_upsample(image)
% image_upsample: Upsamples image by factor of 2 using FFT
% usage: imnew = image_upsample(image)
%   Image dimensions must be powers of 2 

[m n] = size(image);

F = fftshift(fft2(image));
PF = padarray(F, [m/2, n/2]);
imnew = ifft2(ifftshift(PF),'symmetric');
