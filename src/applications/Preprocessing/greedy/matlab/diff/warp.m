function J = warp(I, wx, wy)

[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));
J = interp2(I, mx + wx, my + wy, 'linear', 0);