function [fx fy dfx dfy] = expmap_gateaux_deriv(n, wx, wy, dwx, dwy)
% expmap_gateaux_deriv - compute Gateaux derivative of phi 
% usage: y = expmap(n, wx, wy)
%   n          Number of steps (7-10 is good)
%   wx,wy      Stationary vector field
%   dwx,dwy    Derivative of the vector field w/r/t given variation
% assumptions:
%   - The spacing between voxels is 1

fx = wx / 2^n;
fy = wy / 2^n;

dfx = dwx / 2^n;
dfy = dwy / 2^n;

[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));

for it=1:n
    % Compute psi(t/2, x + psi(t/2,x))
    nfx = fx + interp2(fx, mx + fx, my + fy, '*linear', 0);
    nfy = fy + interp2(fy, mx + fx, my + fy, '*linear', 0);
    
    % Compute Identity plus Jacobian of psi(t/2, *) at x+psi(t/2,x); 
    Jxx = interp2(fx, mx + fx + 0.5, my + fy, '*nearest', 0) - ...
        interp2(fx, mx + fx - 0.5, my + fy, '*nearest', 0);
    Jyx = interp2(fy, mx + fx + 0.5, my + fy, '*nearest', 0) - ...
        interp2(fy, mx + fx - 0.5, my + fy, '*nearest', 0);
    Jxy = interp2(fx, mx + fx, my + fy + 0.5, '*nearest', 0) - ...
        interp2(fx, mx + fx, my + fy - 0.5, '*nearest', 0);
    Jyy = interp2(fy, mx + fx, my + fy + 0.5, '*nearest', 0) - ...
        interp2(fy, mx + fx, my + fy - 0.5, '*nearest', 0);
            
    % [Jxx Jxy] = gradient(fx); [Jyx Jyy] = gradient(fy);
    % Hxx = interp2(Jxx, mx + fx, my + fy, '*linear', 0);
    % Hxy = interp2(Jxy, mx + fx, my + fy, '*linear', 0);
    % Hyx = interp2(Jyx, mx + fx, my + fy, '*linear', 0);
    % Hyy = interp2(Jyy, mx + fx, my + fy, '*linear', 0);
    
    Tx = interp2(dfx, mx + fx, my + fy, '*linear', 0);
    Ty = interp2(dfy, mx + fx, my + fy, '*linear', 0);
    
    % Multiply through to get d(psi(t))
    ndfx = Jxx .* dfx + Jxy .* dfy + Tx + dfx;
    ndfy = Jyx .* dfx + Jyy .* dfy + Ty + dfy;

    % Assign the new variables
    fx = nfx; fy = nfy;
    dfx = ndfx; dfy = ndfy;

    % gridplot(fx,fy,5,5);
    % getframe;
end    
    