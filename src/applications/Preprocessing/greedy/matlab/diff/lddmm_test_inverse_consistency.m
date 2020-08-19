function lddmm_test_inverse_consistency(p)
% lddmm_test_inverse_consistency - test inverse of diffeomorphism
% usage: lddmm_test_inverse_consistency(p)

% Create a random field
[tvx tvy] = lddmm_random_field(p,3);

% Integrate the field
[ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(tvx, tvy, p);

% Get the diffs we care about
f10x = ft0x(:,:,end); f10y = ft0y(:,:,end);
f01x = ft1x(:,:,1);   f01y = ft1y(:,:,1);

% Compute residual
[resx resy] = lddmm_compose_diffeomorphisms(f10x, f10y, f01x, f01y, p);
resabs = sqrt(resx.^2+resy.^2);

% Show grid plot of residual
subplot(1,2,1);
gridplot(resx, resy, 5, 5);

% Print max/RMS statistics
fprintf('Residual Max = %d, RMS = %d\n', ...
    max(resabs(:)), sqrt(mean(mean(resabs.^2))));

% Compute residual
[resx resy] = lddmm_compose_diffeomorphisms(f01x, f01y, f10x, f10y, p);
resabs = sqrt(resx.^2+resy.^2);

% Show grid plot of residual
subplot(1,2,2);
gridplot(resx, resy, 5, 5);

% Print max/RMS statistics
fprintf('Residual Max = %d, RMS = %d\n', ...
    max(resabs(:)), sqrt(mean(mean(resabs.^2))));