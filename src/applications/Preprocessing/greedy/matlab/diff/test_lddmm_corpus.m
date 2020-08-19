
img_c=padarray(double(imread('oval.png')) / 256,[14 14]);
img_o=padarray(double(imread('corpus.png')) / 256,[14 14]);

% Smooth and resample to 128x128
sim_c = imfilter(img_c,fspecial('gaussian',33,1));
sim_o = imfilter(img_o,fspecial('gaussian',33,1));

% Create initial problem 
% [vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.01 / 4 * size(sim_o,1) * size(sim_o,2), 1, 0.008);
[vx vy p] = lddmm_init(sim_o,sim_c, 0.1, 0.001 * size(sim_o,1) * size(sim_o,2), 0.2, 0.004);

% Stick the contour into p'
cont = 0.5 * max(max(sim_o));
p.contour = cont;

% Optimize
% [vxo2 vyo2] = lddmm_optimize_lbgfs(vx,vy,p,60);
[vxo2 vyo2] = lddmm_optimize_basic(vx,vy,p,0.2e-4,200);


