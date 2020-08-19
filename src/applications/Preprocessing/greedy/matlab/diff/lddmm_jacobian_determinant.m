function J = lddmm_jacobian_determinant(fx, fy, p)
% lddmm_jacobian_determinant - compute Jacobian determinant of a map
% usage: 
%   J = lddmm_jacobian_determinant(fx, fy, p)
% params:
%   fx, fy      Displacement fields
%   p           See lddmm_init
[fxx fxy] = gradient(p.mx + fx);
[fyx fyy] = gradient(p.my + fy);
J = fxx .* fyy - fxy .* fyx;
