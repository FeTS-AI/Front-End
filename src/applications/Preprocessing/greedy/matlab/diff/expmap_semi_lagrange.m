function [fx fy] = expmap_semi_lagrange(dt, wx, wy)

% Initialize grid
[m n] = size(wx);
[mx my] = meshgrid(1:n,1:m);

% Set fx, fy at time zero
fx = mx; fy = my;

% Iterate over time steps
for t=dt:dt:1
    
    % Estimate alpha
    ax = zeros(m,n); ay = zeros(m,n);
    for i = 1:5
        nax = dt * interp2(wx,mx-ax/2,my-ay/2,'*linear',0); 
        nay = dt * interp2(wy,mx-ax/2,my-ay/2,'*linear',0);
        ax = nax; ay = nay;
    end
    
    % Interpolate previous diff at alpha
    nfx = interp2(fx, mx-ax, my-ay, '*linear', 0);
    nfy = interp2(fy, mx-ax, my-ay, '*linear', 0);
    fx = nfx; fy = nfy;
end

fx = fx - mx; fy = fy - my;