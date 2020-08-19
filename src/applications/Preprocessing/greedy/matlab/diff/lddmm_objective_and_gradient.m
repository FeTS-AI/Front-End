function [E dedvx dedvy opts] = lddmm_objective_and_gradient(vx, vy, p)
% lddmm_objective_and_gradient - compute LDDMM objective and gradient
% usage: 
%   E = lddmm_objective_and_gradient(vx, vy, p)
%   [E dedvx dedvy] = lddmm_objective_and_gradient(vx, vy, p)
%   [E dedvx dedvy opts] = lddmm_objective_and_gradient(vx, vy, p)
% returns:
%   opts        Specify to return some internal structures

% Compute the diffeomorphisms
[ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p);

% Flow image IO to time 1
J0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,end), ft0y(:,:,end), p);

% Compute the field norm component
e_field = lddmm_vector_field_dot_product(vx, vy, vx, vy, p);

% Compute the image norm component
e_image = sum(sum((J0 - p.I1).^2));

% Compute the total energy 
E = e_field + e_image / p.sigma^2;

%fprintf('Energy components %f, %f\n', e_field, e_image / p.sigma^2);

% Now, compute the gradient if required
if nargout > 1
   
    % The field component of the gradient is trivial
    dedvx = 2 * vx; dedvy = 2 * vy;
    
    % Now the part of the kernel. Get image J_t^0 and J_t^1
    for it = 1:p.nt
        
        % Warp images I0 and I1 to the current time point
        Jt0 = lddmm_warp_scalar_field(p.I0, ft0x(:,:,it), ft0y(:,:,it), p); 
        Jt1 = lddmm_warp_scalar_field(p.I1, ft1x(:,:,it), ft1y(:,:,it), p); 
        
        % Take the gradient of the warped I0 image
        [grad_Jt0_x grad_Jt0_y] = gradient(Jt0);

        % Scale the gradient by the determinant of the (t,1) warp to
        % account for the difference in metric tensor
        detjac_phi_t1 = lddmm_jacobian_determinant(ft1x(:,:,it), ft1y(:,:,it), p);
                
        % Get the right hand side of the PDE
        pde_rhs_x = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_x; 
        pde_rhs_y = detjac_phi_t1 .* (Jt0 - Jt1) .* grad_Jt0_y; 
        
        % Apply the PDE using fft
        pde_soln_x = ifft2(fft2(pde_rhs_x) ./ p.f_kernel_sq,'symmetric');
        pde_soln_y = ifft2(fft2(pde_rhs_y) ./ p.f_kernel_sq,'symmetric');
        
        % Add the solution to the gradient
        dedvx(:,:,it) = dedvx(:,:,it) - 2 * pde_soln_x / p.sigma^2;
        dedvy(:,:,it) = dedvy(:,:,it) - 2 * pde_soln_y / p.sigma^2;        
                
    end
end

% Store the opts array
if nargout > 3
    opts=struct('ft0x',ft0x,'ft0y',ft0y,'ft1x',ft1x,'ft1y',ft1y);
end