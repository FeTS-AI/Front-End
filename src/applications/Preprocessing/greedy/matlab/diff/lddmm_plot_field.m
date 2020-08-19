function lddmm_plot_field(vx, vy)
% lddmm_plot_field(vx, vy) - plot time-varying field

nt = size(vx,3);
for i = 1:nt
    subplot(ceil(sqrt(nt)),ceil(sqrt(nt)),i);
    imagesc(vx(:,:,i));
    axis image;
end