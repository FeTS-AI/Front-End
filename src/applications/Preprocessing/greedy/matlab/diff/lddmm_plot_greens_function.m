function lddmm_plot_greens_function(p)
% lddmm_plot_greens_function - plot green's function of the hilbert space
%   kernel specified by alpha and gamma
% usage:
%   lddmm_plot_greens_function(p)

z = zeros(p.m, p.n);
z(round(p.m / 2), round(p.n / 2)) = 1;
g1 = ifft2(fft2(z) ./ p.f_kernel,'symmetric');
g2 = ifft2(fft2(z) ./ p.f_kernel_sq,'symmetric');
subplot(1,2,1);
imagesc(g1); axis image; title('Green Fn of \gamma I = \alpha \nabla'); colorbar;
subplot(1,2,2);
imagesc(g2); axis image; title('Green Fn of (\gamma I = \alpha \nabla)^2');  colorbar;
