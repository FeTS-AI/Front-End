% Define a random deformation field
rx=padarray(random('norm',0,10,156,156),[50 50]);
ry=padarray(random('norm',0,10,156,156),[50 50]);
wx=20 * imfilter(rx, fspecial('gaussian',80,10));
wy=20 * imfilter(ry, fspecial('gaussian',80,10));
[mx my] = meshgrid(1:size(wx,1),1:size(wx,2));

% Compute forward and inverse transforms
[fx fy] = expmap(7,wx,wy);
[gx gy] = expmap(7,-wx,-wy);

% Compose forward and inverse
dx = warp(gx, fx, fy) + fx;
dy = warp(gy, fx, fy) + fy;

% Look the residuals
fprintf('Inverse Residuals: %f, %f\n', ...
    sqrt(sum(sum(dx.^2)) / length(dx(:))), ...
    sqrt(sum(sum(dy.^2)) / length(dy(:))));

% Look at image registration
img = double(imread('diff_c.png'));
img2 = warp(warp(img, fx, fy), gx, gy);

clf;
gridplot(dx, dy, 5, 5);