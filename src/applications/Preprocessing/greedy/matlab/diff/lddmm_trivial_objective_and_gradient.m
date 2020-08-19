function [E dedvx dedvy opts] = lddmm_trivial_objective_and_gradient(vx, vy, p)
% lddmm_trivial_objective_and_gradient - compute objective and gradient
% of a simple functional E(v) = |phi_10|^2
% usage: 
%   E = lddmm_objective_and_gradient(vx, vy, p)
%   [E dedvx dedvy] = lddmm_objective_and_gradient(vx, vy, p)
%   [E dedvx dedvy opts] = lddmm_objective_and_gradient(vx, vy, p)
% returns:
%   opts        Specify to return some internal structures

% Compute the diffeomorphisms
if nargout > 1
    [ft0x ft0y ft1x ft1y j] = ...
        lddmm_integrate_field_and_jacobian_semi_lagrangian(vx, vy, p);
else
    [ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p);
end

% Compute the total energy 
E = sum(sum(ft0x(:,:,end).^2 + ft0y(:,:,end).^2));

% Now, compute the gradient if required
if nargout > 1
   
    % Initialize at zero
    dedvx = zeros(size(vx)); dedvy = zeros(size(vy));
    
    % Compute the variation
    for it = 1:p.nt
        
        % Compute the Jacobian of ft0
        fxx = j.Dft0_xx(:,:,it);
        fxy = j.Dft0_xy(:,:,it);
        fyx = j.Dft0_yx(:,:,it);
        fyy = j.Dft0_yy(:,:,it);
        
        % Compute the Jacobian determinant of ft1
        g = j.Dft1_xx(:,:,it) .* j.Dft1_yy(:,:,it) - ...
            j.Dft1_xy(:,:,it) .* j.Dft1_yx(:,:,it);
        
        % Compute the result
        dedvx(:,:,it) = -2 * (g .* (fxx .* ft0x(:,:,it) + fxy .* ft0y(:,:,it)));
        dedvy(:,:,it) = -2 * (g .* (fyx .* ft0x(:,:,it) + fyy .* ft0y(:,:,it)));
                
    end
end

% Store the opts array
if nargout > 3
    opts=struct('ft0x',ft0x,'ft0y',ft0y,'ft1x',ft1x,'ft1y',ft1y);
end