function [e_total dwx dwy] = diff_objective(wx, wy, p)

I0 = p.I0;
I1 = p.I1;
mx = p.mx;
my = p.my;
A = p.A;
A2 = p.A2;
sigma = p.sigma;

% Compute \phi^-1 - Id 
[fx0 fy0] = expmap(7, -wx, -wy);

% Compute I_0 \circ \phi^-1 
J0 = interp2(I0, mx + fx0, my + fy0, '*linear', 0);

% Compute \nabla(I_0 \circ \phi^{-1})
[I0dx I0dy] = gradient(I0);
[fxdx fxdy] = gradient(fx0);
[fydx fydy] = gradient(fy0);
Gx0 = I0dx .* (ones(size(I0)) + fxdx) + I0dy .* fydx;
Gy0 = I0dx .* fxdy + I0dy .* (ones(size(I0)) + fydy);

%[Gx0 Gy0] = gradient(J0);

% Compute the term on the right hand of L L-adj in (12)
dIx = (J0 - I1) .* Gx0;
dIy = (J0 - I1) .* Gy0;

% Compute the update step
dwx = -2 * dIx;
dwy = -2 * dIy;

% Compute the image match energy
e_image = sum(sum((J0 - I1).^2));

% Print the energy
e_total = e_image;

e_total = sum(sum(fx0)) + sum(sum(fy0));
dwx = -2 * wx;
dwy = -2 * wy;
