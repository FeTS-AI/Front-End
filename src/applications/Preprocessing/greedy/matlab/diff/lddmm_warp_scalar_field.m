function J = lddmm_warp_scalar_field(I, fx, fy, p)
% lddmm_warp_scalar_field - warp scalar field
% usage: 
%   J = lddmm_warp_scalar_field(I, fx, fy, p)
% params:
%   I       : scalar field
%   fx, fy  : displacement field (2D arrays)
%   p       : see lddmm_init

J = interp2(I, p.mx + fx, p.my + fy, '*linear', 0);