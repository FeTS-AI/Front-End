function [dx dy vxo2 vyo2] = lddmm_basic(img_fix, img_mov, nb_iter, vx1, vy1)% Load the images
% gen_roof_image
% img_c = I1 / 200;
% img_o = I2 / 200;
% 
sim_c=img_mov; 
sim_o=img_fix;
% img_c=double(imread('diff_circ2.png')) / 256;
% img_o=double(imread('diff_circ.png')) / 256;
% Smooth and resample to 128x128
% sim_c = imresize(imfilter(img_c,fspecial('gaussian',32,4)),[128 128]);
% sim_o = imresize(imfilter(img_o,fspecial('gaussian',32,4)),[128 128]);

% sim_c = imresize(img_c, [128 128]);
% sim_o = imresize(img_o, [128 128]);

% Create initial problem 
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * size(sim_o,1) * size(sim_o, 2), 1, 0.008);
% [vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 * 64 * 64, 1, 0.001);
% [vx vy p] = lddmm_init(sim_o,sim_c, 0.2, 0.01, 1, 1);


if exist('vx1', 'var')
    vx = vx1;
end;

if exist('vy1', 'var')
    vy = vy1;
end;

% Optimize
% [vxo2 vyo2] = lddmm_optimize_lbgfs(vx,vy,p,60);
[vxo2 vyo2 q] = lddmm_optimize_basic(vx,vy,p,1e-2,nb_iter);


dx = -1 * q.ft0x(:,:,end);
dy = -1 * q.ft0y(:,:,end);