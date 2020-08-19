function [dwx dwy] = elupdate(I0, I1, wx, wy, gamma, alpha, sigma)

[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));

[fx0 fy0] = expmap(7, -wx, -wy);
[fx1 fy1] = expmap(7, wx, wy);

J0 = interp2(I0, mx + fx0, my + fy0, '*linear', 0);
J1 = interp2(I1, mx + fx1, my + fy1, '*linear', 0);

[Gx0 Gy0] = gradient(J0);
[Gx1 Gy1] = gradient(J1);

dIx = (J0 - I1) .* Gx0 + (J1 - I0) .* Gx1;
dIy = (J0 - I1) .* Gy0 + (J1 - I0) .* Gy1;

Zx = fft2(dIx);
Zy = fft2(dIy);

A = gamma + 2 * alpha * (2 - cos(2 * pi * (mx-1) / size(wx,1)) - cos(2 * pi * (my-1) / size(wx,2)));
A2 = A.*A;
LLIx = real(ifft2(Zx ./ A2));
LLIy = real(ifft2(Zy ./ A2));

dwx = 2 * (wx - LLIx / sigma^2);
dwy = 2 * (wy - LLIy / sigma^2);
