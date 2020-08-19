function [fx fy] = expmap(n, wx, wy)
% expmap - compute exponential map of 2d vector field
% usage: y = expmap(n, wx, wy)
%   n        Number of steps (7-10 is good)
%   wx,wy    Stationary vector field
%   fx,fy    Displacement field at time t 
% assumptions:
%   - The spacing between voxels is 1

fx = wx / 2^n;
fy = wy / 2^n;

[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));

for it=1:n
    ndx = interp2(fx, mx + fx, my + fy, '*linear', 0);
    ndy = interp2(fy, mx + fx, my + fy, '*linear', 0);
    fx = fx + ndx;
    fy = fy + ndy;

    %gridplot(fx,fy,5,5);
    %drawnow
end    
    