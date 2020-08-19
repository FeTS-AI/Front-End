function [ft0x ft0y ft1x ft1y j] = ...
    lddmm_integrate_field_and_jacobian_semi_lagrangian(vx, vy, p)
% lddmm_integrate_field_semi_lagrangian - integrate a diffeomorphism from 
%   a time-dependent vector field
% usage:
%   [ft0x ft0y ft1x ft1y j ] = 
%        lddmm_integrate_field_and_jacobian_semi_lagrangian(vx, vy, p)
% parameters:
%   vx, vy:     Time-varying vector field V
%   p           Parameter structure
% returns:
%   ft0x, ft0y:   Transform phi_{t_j,0} from time t to time 0 (3d array) 
%   ft1x, ft1y:   Transform phi_(t_j,1) from time t to time 1 (3d array)
%   j             Struct where the Jacobian will be placed

% Initialize all transforms to zero
ft0x = zeros(p.m, p.n, p.nt);
ft0y = zeros(p.m, p.n, p.nt);
ft1x = zeros(p.m, p.n, p.nt);
ft1y = zeros(p.m, p.n, p.nt);

% Initialize the Jacobians
j.Dft0_xx = ones(p.m, p.n, p.nt);
j.Dft0_xy = zeros(p.m, p.n, p.nt);
j.Dft0_yx = zeros(p.m, p.n, p.nt);
j.Dft0_yy = ones(p.m, p.n, p.nt);

j.Dft1_xx = ones(p.m, p.n, p.nt);
j.Dft1_xy = zeros(p.m, p.n, p.nt);
j.Dft1_yx = zeros(p.m, p.n, p.nt);
j.Dft1_yy = ones(p.m, p.n, p.nt);

% Allocate an array of alphas (not to compute alphas twice)
ax = zeros(p.m, p.n, p.nt);
ay = zeros(p.m, p.n, p.nt);
Da_xx = zeros(p.m, p.n, p.nt);
Da_xy = zeros(p.m, p.n, p.nt);
Da_yx = zeros(p.m, p.n, p.nt);
Da_yy = zeros(p.m, p.n, p.nt);

for it = 1:p.nt
    
    [Dv_xx Dv_xy] = gradient(vx(:,:,it));
    [Dv_yx Dv_yy] = gradient(vy(:,:,it));
    
    for i = 1:5
        nax = p.dt * interp2(vx(:,:,it), p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2,'*linear',0); 
        nay = p.dt * interp2(vy(:,:,it), p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2,'*linear',0);
        
        txx = p.dt * interp2(Dv_xx, p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2, '*linear',0);
        txy = p.dt * interp2(Dv_xy, p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2, '*linear',0);
        tyx = p.dt * interp2(Dv_yx, p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2, '*linear',0);
        tyy = p.dt * interp2(Dv_yy, p.mx-ax(:,:,it)/2, p.my-ay(:,:,it)/2, '*linear',0);
        
        qxx = txx - 0.5 * (Da_xx(:,:,it) .* txx - Da_yx(:,:,it) .* txy);
        qxy = txy - 0.5 * (Da_xy(:,:,it) .* txx - Da_yy(:,:,it) .* txy);
        qyx = tyx - 0.5 * (Da_xx(:,:,it) .* tyx - Da_yx(:,:,it) .* tyy);
        qyy = tyy - 0.5 * (Da_xy(:,:,it) .* tyx - Da_yy(:,:,it) .* tyy);
        
        ax(:,:,it) = nax; ay(:,:,it) = nay;
        
        Da_xx(:,:,it) = qxx;
        Da_xy(:,:,it) = qxy;
        Da_yx(:,:,it) = qyx;
        Da_yy(:,:,it) = qyy;
    end  
end

% Iterate over time steps
for it = 2:p.nt
    % Interpolate previous diff at alpha
    ft0x(:,:,it) = interp2(ft0x(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear', 0) - ax(:,:,it);
    ft0y(:,:,it) = interp2(ft0y(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear', 0) - ay(:,:,it);
    
    % Compute the Jacobians
    txx = interp2(j.Dft0_xx(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear',1);
    txy = interp2(j.Dft0_xy(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear',0);
    tyx = interp2(j.Dft0_yx(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear',0);
    tyy = interp2(j.Dft0_yy(:,:,it-1), p.mx-ax(:,:,it), p.my-ay(:,:,it), '*linear',1);
    
    j.Dft0_xx(:,:,it) = txx - Da_xx(:,:,it) .* txx - Da_yx(:,:,it) .* txy;
    j.Dft0_xy(:,:,it) = txy - Da_xy(:,:,it) .* txx - Da_yy(:,:,it) .* txy;
    j.Dft0_yx(:,:,it) = tyx - Da_xx(:,:,it) .* tyx - Da_yx(:,:,it) .* tyy;
    j.Dft0_yy(:,:,it) = tyy - Da_xy(:,:,it) .* tyx - Da_yy(:,:,it) .* tyy;
    
end

% Iterate over time steps
for it = p.nt-1:-1:1
    % Interpolate previous diff at alpha
    ft1x(:,:,it) = interp2(ft1x(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear', 0) + ax(:,:,it);
    ft1y(:,:,it) = interp2(ft1y(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear', 0) + ay(:,:,it);
    
    % Compute the Jacobians
    txx = interp2(j.Dft1_xx(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear',1);
    txy = interp2(j.Dft1_xy(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear',0);
    tyx = interp2(j.Dft1_yx(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear',0);
    tyy = interp2(j.Dft1_yy(:,:,it+1), p.mx+ax(:,:,it), p.my+ay(:,:,it), '*linear',1);
    
    j.Dft1_xx(:,:,it) = txx + Da_xx(:,:,it) .* txx + Da_yx(:,:,it) .* txy;
    j.Dft1_xy(:,:,it) = txy + Da_xy(:,:,it) .* txx + Da_yy(:,:,it) .* txy;
    j.Dft1_yx(:,:,it) = tyx + Da_xx(:,:,it) .* tyx + Da_yx(:,:,it) .* tyy;
    j.Dft1_yy(:,:,it) = tyy + Da_xy(:,:,it) .* tyx + Da_yy(:,:,it) .* tyy;
end

