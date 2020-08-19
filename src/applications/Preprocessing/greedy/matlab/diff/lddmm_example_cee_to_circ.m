% Load the images
img_c=double(imread('diff_circ2.png')) / 256;
img_o=double(imread('diff_circ.png')) / 256;

% Smooth and resample to 128x128
sim_c = imresize(imfilter(img_c,fspecial('gaussian',32,4)),[128 128]);
sim_o = imresize(imfilter(img_o,fspecial('gaussian',32,4)),[128 128]);

% Create initial problem 
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * 64 * 64, 1, 0.008);

% Optimize
[vxo2 vyo2] = lddmm_optimize_lbgfs(vx,vy,p,60);