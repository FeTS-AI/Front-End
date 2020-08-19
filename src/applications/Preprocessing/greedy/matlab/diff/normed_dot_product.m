function f = normed_dot_product(ax, ay, bx, by, p)

Lax = ifft2(p.A2 .* fft2(ax),'symmetric');
Lay = ifft2(p.A2 .* fft2(ay),'symmetric');
dp = Lax .* bx + Lay .* by;
f = sum(dp(:));