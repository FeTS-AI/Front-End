function [nvx nvy] = lddmm_reparameterize_vector_field(vx, vy, p)
% lddmm_reparameterize_vector_field - reparamterize time-dep vec field
% usage:
%   [nvx nvy] = lddmm_reparameterize_vector_field(vx, vy, p)

% Compute norm of V at each time point
for i = 1:p.nt    
    nv(i) = sqrt(lddmm_hilbert_dot_product(vx(:,:,i), vy(:,:,i), vx(:,:,i), vy(:,:,i), p));
end

plot(nv);

% Integrate norm of V from 0 to 1, normalize to 1
len = trapz(0:p.dt:1,nv);
s = cumtrapz(0:p.dt:1,nv) / len;

% Compute the inverse of s(t), h
h = zeros(size(s));
h(end) = 1;
for i = 2:p.nt-1
    t = p.dt * (i-1);
    h(i) = fminbnd(@(x)((t - interp1(0:p.dt:1,s,[x])))^2,0,1);
end

% Resample vector fields at h(t)
nvx = shiftdim(interp1(0:p.dt:1, shiftdim(vx,2), h),1);
nvy = shiftdim(interp1(0:p.dt:1, shiftdim(vy,2), h),1);

% Compute h'(t)
hprime = interp1(0:p.dt:1, nv/len, h) .^ -1;

% Scale vector fields by h'
for i = 1:p.nt
    nvx(:,:,i) = hprime(i) * nvx(:,:,i);
    nvy(:,:,i) = hprime(i) * nvy(:,:,i);
end

for i = 1:p.nt    
    nnv(i) = sqrt(lddmm_hilbert_dot_product(nvx(:,:,i), nvy(:,:,i), nvx(:,:,i), nvy(:,:,i), p));
end

plot(nv);
hold on;
plot(nnv,'r');
hold off;
axis tight;