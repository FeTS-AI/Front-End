% Load the two images
img_c=double(imread('brain_fixed.png')) / 256;
img_o=double(imread('brain_moving.png')) / 256;

% Smooth and resample to 128x128
sim_c = imresize(imfilter(img_c,fspecial('gaussian',13,1.2)),[128 128]);
sim_o = imresize(imfilter(img_o,fspecial('gaussian',13,1.2)),[128 128]);

% Extract the contour from the fixed image
[idx ctr] = kmeans(sim_o(sim_o > 0), 3);
cont = conv(sort(ctr), [0.5 0.5], 'valid');

% Create initial problem 
% [vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 / 4 * size(sim_o,1) * size(sim_o,2), 1, 0.008);
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * size(sim_o,1) * size(sim_o,2), 1, 0.0001);

% Stick the contour into p
p.contour = cont;

% Optimize
% [vxo2 vyo2] = lddmm_optimize_lbgfs(vx,vy,p,60);
[vxo2 vyo2] = lddmm_optimize_basic(vx,vy,p,1e-6,50);


