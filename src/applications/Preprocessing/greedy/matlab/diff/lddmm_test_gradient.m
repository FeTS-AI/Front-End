function lddmm_test_gradient(vx,vy,p)
% lddmm_test_gradient - test gradient computation
% usage:
%   lddmm_test_gradient(vx,vy,p)
% parameters:
%   vx,vy   Vector field to evaluate at
%   p       Generated using lddmm_init

% Evaluate at current point
[E dedvx dedvy] = lddmm_objective_and_gradient(vx, vy, p);

% Compute the Euclidean gradient
gex = dedvx;
gey = dedvy;

for i = 1:p.nt
    gex(:,:,i) = ifft2(fft2(dedvx(:,:,i)) .* p.f_kernel_sq, 'symmetric');
    gey(:,:,i) = ifft2(fft2(dedvy(:,:,i)) .* p.f_kernel_sq, 'symmetric');
end 

%varx = gex; vary = gey;
varx = dedvx; vary = dedvy;

% Loop over random variations
for i = 1:10
    
    % gateaux_analytic = lddmm_vector_field_dot_product(dedvx, dedvy, varx, vary, p);
    gateaux_analytic = sum(sum(sum(gex .* varx + gey.* vary))) / p.nt;
    
    for expeps = -6:-6
        
        eps = 10^expeps;
    
        E1 = lddmm_objective_and_gradient(vx - eps * varx, vy - eps * vary, p);
        E2 = lddmm_objective_and_gradient(vx + eps * varx, vy + eps * vary, p);

        gateaux_numeric = (E2 - E1) / (2 * eps);

        fprintf('Iter: %4i    Eps: %4.2e    Analytic: %12d     Numeric: %12d\n', ...
            i, eps, gateaux_analytic, gateaux_numeric);
        
    end
    
    [varx vary] = lddmm_random_field(p);
    for i = 1:p.nt
        ker = fspecial('gaussian', [p.m p.n], 20); 
        ker = ker / max(ker(:));
        varx(:,:,i) = 20 * varx(:,:,i) .* ker;
        vary(:,:,i) = 20 * vary(:,:,i) .* ker;
        imagesc(varx(:,:,i)); axis image; colorbar; caxis([-1 1]); getframe;
    end
end