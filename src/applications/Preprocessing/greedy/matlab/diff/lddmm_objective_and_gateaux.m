function [E dedvx dedvy] = lddmm_objective_and_gradient(vx, vy, p)
% lddmm_objective_and_gradient - compute LDDMM objective and gradient
% usage: 
%   E = lddmm_objective_and_gradient(vx, vy, p)
%   [E dedvx dedvy] = lddmm_objective_and_gradient(vx, vy, p)

% Delete this
fun1=@(x,y)(sin(2*pi*x/128) .* cos(2 * pi * y/64));
fun2=@(x,y)(sin(2*pi*x/64) .* cos(2 * pi * y/128) + sin(2*pi*x/32) .* cos(2 * pi * y/32));
gfun1=@(x,y)((2*pi*cos(2*pi*x/128)/128).*cos(2*pi*y/64));
gfun2=@(x,y)(-sin(2*pi*x/128).*(2*pi*sin(2*pi*y/64)));

% Compute the diffeomorphisms
[ffx ffy fix fiy] = lddmm_integrate_field_semi_lagrangian(vx, vy, p);

% Apply inverse warp to image I0

J0 = fun1(p.mx + fix(:,:,end), p.my + fiy(:,:,end));
% J0 = interp2(p.I0, p.mx + fix(:,:,end), p.my + fiy(:,:,end), '*linear', 0);

% Compute the field norm component
e_field = lddmm_vector_field_dot_product(vx, vy, vx, vy, p);

% Compute the image norm component
e_image = sum(sum((J0 - p.I1).^2));

% Compute the total energy 
E = e_field + e_image / p.sigma^2;

% Now, compute the gradient if required
if nargout > 1
   
    % The field component of the gradient is trivial
    dedvx = 2 * vx; dedvy = 2 * vy;
    
    % Now the part of the kernel. Get image J_t^0 and J_t^1
    for it = 1:p.nt
        
        % Compute the right hand side of the PDE 
        % Jt0 = interp2(p.I0, p.mx + fix(:,:,it), p.my + fiy(:,:,it), '*linear', 0);
        % Jt1 = interp2(p.I1, p.mx + ffx(:,:,it), p.my + ffy(:,:,it), '*linear', 0);
        
        Jt0 = fun1(p.mx + fix(:,:,it), p.my + fiy(:,:,it));
        Jt1 = fun2(p.mx + ffx(:,:,it), p.my + ffy(:,:,it));
        
        %[grad_Jt0_x grad_Jt0_y] = gradient(Jt0);
        grad_Jt0_x = gfun1(p.mx + fix(:,:,it), p.my + fiy(:,:,it));
        grad_Jt0_y = gfun2(p.mx + fix(:,:,it), p.my + fiy(:,:,it));
        
        [dffx_dx dffx_dy] = gradient(p.mx + ffx(:,:,it));
        [dffy_dx dffy_dy] = gradient(p.my + ffy(:,:,it));
        detjac_ff = dffx_dx .* dffy_dy - dffx_dy .* dffy_dx;
        pde_rhs_x = detjac_ff .* (Jt0 - Jt1) .* grad_Jt0_x; 
        pde_rhs_y = detjac_ff .* (Jt0 - Jt1) .* grad_Jt0_y; 
        
        % Apply the PDE using fft
        pde_soln_x = ifft2(fft2(pde_rhs_x) ./ p.f_kernel_sq,'symmetric');
        pde_soln_y = ifft2(fft2(pde_rhs_y) ./ p.f_kernel_sq,'symmetric');
        
        % Add the solution to the gradient
        dedvx(:,:,it) = dedvx(:,:,it) - 2 * pde_soln_x / p.sigma^2;
        dedvy(:,:,it) = dedvy(:,:,it) - 2 * pde_soln_y / p.sigma^2;        
        
    end
end
