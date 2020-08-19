function [ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p)
% lddmm_integrate_field_semi_lagrangian - integrate a diffeomorphism from 
%   a time-dependent vector field
% usage:
%   [ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p)
% parameters:
%   vx, vy:     Time-varying vector field V
%   p           Parameter structure
% returns:
%   ft0x, ft0y:   Transform phi_{t_j,0} from time t to time 0 (3d array) 
%   ft1x, ft1y:   Transform phi_(t_j,1) from time t to time 1 (3d array)

% Initialize all transforms to zero
ft0x = zeros(p.m, p.n, p.nt);
ft0y = zeros(p.m, p.n, p.nt);
ft1x = zeros(p.m, p.n, p.nt);
ft1y = zeros(p.m, p.n, p.nt);

% Allocate an array of alphas (not to compute alphas twice)
ax = zeros(p.m, p.n, p.nt);
ay = zeros(p.m, p.n, p.nt);
for it = 1:p.nt
    for i = 1:5
        nax = p.dt * interp2(vx(:,:,it), p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2,'*linear',0); 
        nay = p.dt * interp2(vy(:,:,it), p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2,'*linear',0);
        ax(:,:,it) = nax; ay(:,:,it) = nay;
    end  
end



% Iterate over time steps
for it = 2:p.nt
    % Interpolate previous diff at alpha
    ft0x(:,:,it) = interp2(ft0x(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear', 0) - ax(:,:,it);
    ft0y(:,:,it) = interp2(ft0y(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear', 0) - ay(:,:,it);
end

% Iterate over time steps
for it = p.nt-1:-1:1
    % Interpolate previous diff at alpha
    ft1x(:,:,it) = interp2(ft1x(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear', 0) + ax(:,:,it);
    ft1y(:,:,it) = interp2(ft1y(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear', 0) + ay(:,:,it);
end

