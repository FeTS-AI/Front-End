function I = mkgridimage(nx, ny, ix, iy)

I = zeros(nx,ny);
I(1:ix:nx,:) = 1;
I(:,1:iy:ny) = 1;