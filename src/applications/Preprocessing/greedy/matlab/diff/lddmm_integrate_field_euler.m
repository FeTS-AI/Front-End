function [f01x f01y f10x f10y] = lddmm_integrate_field_euler(vx, vy, p)

f01x = zeros(p.m, p.n);
f01y = zeros(p.m, p.n);

f10x = zeros(p.m, p.n);
f10y = zeros(p.m, p.n);

for i = 1:p.nt
    nf01x = f01x + p.dt * interp2(vx(:,:,i), p.mx + f01x, p.my + f01y, 'linear', 0);
    nf01y = f01y + p.dt * interp2(vy(:,:,i), p.mx + f01x, p.my + f01y, 'linear', 0);
    f01x = nf01x;
    f01y = nf01y;

    nf10x = f10x - p.dt * interp2(vx(:,:,1+end-i), p.mx + f10x, p.my + f10y, 'linear', 0);
    nf10y = f10y - p.dt * interp2(vy(:,:,1+end-i), p.mx + f10x, p.my + f10y, 'linear', 0);
    f10x = nf10x;
    f10y = nf10y;

end
    

