function lddmm_plot_problem(vx, vy, p, opts)
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

% Flow image IO to time 1
J0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,end), ft0y(:,:,end), p);
J1 = lddmm_warp_scalar_field(p.I1, ft1x(:,:,1), ft1y(:,:,1), p);

% Plot images
subplot(3,3,1);
imagesc(p.I0); axis image; 
title('I_0');

subplot(3,3,2);
imagesc(J1); axis image; 
title('I_1 \circ \phi_{0,1} (I_1 warped to I_0)');

subplot(3,3,3);
imagesc(J1 - p.I0); caxis([-1 1]);  axis image;  
title('I_1 \circ \phi_{0,1} - I_0 (Not the objective)');

subplot(3,3,4);
imagesc(p.I1); axis image; 
title('I_1');

subplot(3,3,5);
imagesc(J0); axis image; 
title('I_0 \circ \phi_{1,0} (I_0 warped to I_1)');

subplot(3,3,6);
imagesc(J0 - p.I1); caxis([-1 1]); axis image; 
title('I_0 \circ \phi_{1,0} - I_1 (Objective)');

% Plot grids
subplot(3,3,7);
gridplot(ft0x(:,:,end), ft0y(:,:,end),5,5);
title('\phi_{1,0}');

subplot(3,3,8);
gridplot(ft1x(:,:,1), ft1y(:,:,1),5,5);
title('\phi_{0,1}');

% Plot W
subplot(3,3,9);
imagesc(sqrt(mean(vx,3).^2 + mean(vy,3).^2)); axis image; colormap default; caxis([-16 16]);
title('|v|');