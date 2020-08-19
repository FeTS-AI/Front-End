function lddmm_plot_problem_simple(vx, vy, p, opts)
% lddmm_plot_problem - plot state of LDDMM problem
% usage: 
%   lddmm_plot_problem(vx, vy, p, opts)
% parameters:
%   opts        Optional, output of lddmm_objective_and_gradient
%               If specified, plotting takes much less time

% Compute the diffeomorphisms
if nargin < 4
    [ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p);
else
    ft0x = opts.ft0x; ft0y = opts.ft0y;
    ft1x = opts.ft1x; ft1y = opts.ft1y;
end

% Plot the fixed image
subplot(1,2,2);
imagesc(p.I1); colormap gray; axis image;
hold on; 
gridplot(ft0x(:,:,end), ft0y(:,:,end),5,5);
axis([0 size(p.I1,1) 0 size(p.I1,2)]);
title('Fixed Image (I1)');

% Flow image IO to time 1
J0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,end), ft0y(:,:,end), p);
J1 = lddmm_warp_scalar_field(p.I1, ft1x(:,:,1), ft1y(:,:,1), p);

% Plot the moving image
subplot(1,2,1);
imagesc(J0); colormap gray; axis image;
hold on; 
title('Warped Moving Image (I0)');

if isfield(p,'contour')
    hold on; contour(p.I1, p.contour, 'g'); hold off;
end