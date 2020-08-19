function [E dE] = diff_objective_gateux_deriv(wx, wy, hx, hy, p)

I0 = p.I0;
I1 = p.I1;
mx = p.mx;
my = p.my;
A = p.A;
A2 = p.A2;
sigma = p.sigma;

% Compute \phi^-1 - Id 
% Also take Gateaux deriv of the map 
[fx0 fy0 dfx0 dfy0] = expmap_gateaux_deriv(7, -wx, -wy, -hx, -hy);
[fx1 fy1 dfx1 dfy1] = expmap_gateaux_deriv(7, wx, wy, hx, hy);

% Compute I_0 \circ \phi^-1 
J0 = interp2(I0, mx + fx0, my + fy0, '*linear', 0);
Gx0 = ...
    interp2(I0, mx + fx0 + 0.5, my + fy0, '*nearest', 0) - ...
    interp2(I0, mx + fx0 - 0.5, my + fy0, '*nearest', 0);
Gy0 = ...
    interp2(I0, mx + fx0, my + fy0 + 0.5, '*nearest', 0) - ...
    interp2(I0, mx + fx0, my + fy0 - 0.5, '*nearest', 0);
    
% Compute I_1 \circ \phi 
J1 = interp2(I1, mx + fx1, my + fy1, '*linear', 0);
Gx1 = ...
    interp2(I1, mx + fx1 + 0.5, my + fy1, '*nearest', 0) - ...
    interp2(I1, mx + fx1 - 0.5, my + fy1, '*nearest', 0);
Gy1 = ...
    interp2(I1, mx + fx1, my + fy1 + 0.5, '*nearest', 0) - ...
    interp2(I1, mx + fx1, my + fy1 - 0.5, '*nearest', 0);
    
% Compute \nabla(I_0 \circ \phi^{-1})
T0 = Gx0 .* dfx0 + Gy0 .* dfy0;
T1 = Gx1 .* dfx1 + Gy1 .* dfy1;

% Compute the regularization energy
e_field = normed_dot_product(wx, wy, wx, wy, p);
de_field = 2 * normed_dot_product(wx, wy, hx, hy, p);

e_image = (sum(sum((J0 - I1).^2)) + sum(sum((J1 - I0).^2))) / p.sigma^2;
de_image = 2 * (sum(sum((J0 - I1) .* T0)) + sum(sum((J1 - I0) .* T1))) / p.sigma^2;

% Compute the term on the right hand of L L-adj in (12)
dE = de_field + de_image;
E = e_field + e_image;

