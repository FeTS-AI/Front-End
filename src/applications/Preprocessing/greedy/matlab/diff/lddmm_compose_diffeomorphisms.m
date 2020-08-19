function [hx hy] = lddmm_compose_diffeomorphisms(fx, fy, gx, gy, p)
% lddmm_compose_diffeomorphisms 
% usage: 
%   [hx hy] = lddmm_compose_diffeomorphisms(fx, fy, gx, gy)
% notes:
%   all inputs and outputs are displacement fields corresponding
%   to diffeomorphisms F, G, H. The computation is H = F \circ G
% params:
%   fx, fy  : displacement field for F(2D arrays)
%   gx, gy  : displacement field for G(2D arrays)
%   p       : see lddmm_init
hx = gx + interp2(fx, p.mx + gx, p.my + gy, '*linear', 0);
hy = gy + interp2(fy, p.mx + gx, p.my + gy, '*linear', 0);