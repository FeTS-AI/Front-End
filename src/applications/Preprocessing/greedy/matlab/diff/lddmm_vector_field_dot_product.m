function f = lddmm_vector_field_dot_product(ax, ay, bx, by, p)
% lddmm_vector_field_dot_product - computes <A,B>_V, where A and B are
% time-dependent vector fields
% usage:
%   f = lddmm_vector_field_dot_product(ax, ay, bx, by, p)
% parameters:
%   ax,ay   First vector field
%   bx,by   Second vector field
%   p       Created by lddmm_init

f = 0;

for i = 1:p.nt
    f = f + lddmm_hilbert_dot_product(ax(:,:,i), ay(:,:,i), bx(:,:,i), by(:,:,i), p);
end

f = f / p.nt;