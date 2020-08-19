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

% Compute \phi - Id
[fx1 fy1] = expmap(7, wx, wy);

% Compute I_0 \circ \phi^-1 
J0 = interp2(I0, mx + fx0, my + fy0, '*linear', 0);

% Compute the \nabla (I_0 \circ \phi^-1)
Gx0 = ...
    interp2(I0, mx + fx0 + 0.5, my + fy0, '*nearest', 0) - ...
    interp2(I0, mx + fx0 - 0.5, my + fy0, '*nearest', 0);
Gy0 = ...
    interp2(I0, mx + fx0, my + fy0 + 0.5, '*nearest', 0) - ...
    interp2(I0, mx + fx0, my + fy0 - 0.5, '*nearest', 0);

% Compute I_1 \circ \phi 
J1 = interp2(I1, mx - fx1, my - fy1, '*linear', 0);

% Compute the \nabla (I_1 \circ \phi)
Gx1 = ...
    interp2(I1, mx + fx1 + 0.5, my + fy1, '*nearest', 0) - ...
    interp2(I1, mx + fx1 - 0.5, my + fy1, '*nearest', 0);
Gy1 = ...
    interp2(I1, mx + fx1, my + fy1 + 0.5, '*nearest', 0) - ...
    interp2(I1, mx + fx1, my + fy1 - 0.5, '*nearest', 0);

% Compute the component of the gradient corresponding to image forces
dIx = - (J0 - I1) .* Gx0 + (J1 - I0) .* Gx1;
dIy = - (J0 - I1) .* Gy0 + (J1 - I0) .* Gy1;
%dIx = - (J0 - I1) .* Gx0;
%dIy = - (J0 - I1) .* Gy0;

% Apply (LL^{adj})^{-1} to the above
LLIx = ifft2(fft2(dIx) ./ A2, 'symmetric');
LLIy = ifft2(fft2(dIy) ./ A2, 'symmetric');

% Compute the update step
dwx = 2 * (wx + LLIx / sigma^2);
dwy = 2 * (wy + LLIy / sigma^2);

% Compute the image match energy
e_image = (sum(sum((J0 - I1).^2)) + sum(sum((J1 - I0).^2))) / sigma.^2;
%e_image = sum(sum((J0 - I1).^2)) / sigma.^2;

% Compute the regularization energy
e_field = normed_dot_product(wx, wy, wx, wy, p);

% Print the energy
e_total = e_field + e_image;