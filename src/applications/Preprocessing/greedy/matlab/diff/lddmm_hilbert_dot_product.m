function f = lddmm_hilbert_dot_product(ax, ay, bx, by, p)
% lddmm_hilbert_dot_product - compute <a,b>_V = <La,Lb>_{L^2}
% usage:
%   f = lddmm_hilbert_dot_product(ax, ay, bx, by, p)

Lax = ifft2(p.f_kernel_sq .* fft2(ax),'symmetric');
Lay = ifft2(p.f_kernel_sq .* fft2(ay),'symmetric');
dp = Lax .* bx + Lay .* by;
f = sum(dp(:));