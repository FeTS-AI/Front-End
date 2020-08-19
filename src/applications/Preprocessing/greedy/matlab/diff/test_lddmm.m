
img_c=double(imread('circc_fixed.png')) / 256;
img_o=double(imread('circc_moving.png')) / 256;

% Smooth and resample to 128x128
sim_c = imresize(imfilter(img_c,fspecial('gaussian',32,4)),[128 128]);
sim_o = imresize(imfilter(img_o,fspecial('gaussian',32,4)),[128 128]);

% Create initial problem 
% [vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 / 4 * size(sim_o,1) * size(sim_o,2), 1, 0.008);
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * size(sim_o,1) * size(sim_o,2), 1, 0.008);

% Stick the contour into p'
cont = 0.5 * max(max(sim_o));
p.contour = cont;

% Optimize
% [vxo2 vyo2] = lddmm_optimize_lbgfs(vx,vy,p,60);
[vxo2 vyo2] = lddmm_optimize_basic(vx,vy,p,1e-2,200);


