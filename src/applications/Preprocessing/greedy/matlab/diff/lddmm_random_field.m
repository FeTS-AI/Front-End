function [vx vy] = lddmm_random_field(p, n_frames)
% lddmm_random_field - generate random smooth field for LDDMM testing
% usage: [vx vy] = lddmm_random_field(p, n_frames)
% parameters:
%     p         : struct created by lddmm_init
%     n_frames  : number of random 'frames' that are time-interpolated
%                 (value of 1 gives a constant field)

if(nargin < 2)
    n_frames = 3;
end

% Create some noise in 3D
%{
SX=zeros(p.m, p.n, n_frames);
SY=zeros(p.m, p.n, n_frames);

for i = 1:n_frames
    SX(:,:,i) = 20 * imfilter(...
       padarray(random('norm',0,8,[p.m-32,p.n-32]),[16 16]),...
       fspecial('gaussian',80,10));
    SY(:,:,i) = 20 * imfilter(...
       padarray(random('norm',0,8,[p.m-32,p.n-32]),[16 16]),...
       fspecial('gaussian',80,10));
end
%}

SX=zeros(p.m, p.n, n_frames);
SY=zeros(p.m, p.n, n_frames);

for i = 1:n_frames
    for d = 1:2
        K = complex(randn(17), randn(17)) .* fspecial('gaussian', [17 17], 1);
        F = zeros(p.m, p.n);
        cm = floor(p.m/2); cn = floor(p.n/2);
        F(cm-8:cm+8,cn-8:cn+8) = K;
        V = ifft2(fftshift(F), 'symmetric');
        if d==1
            SX(:,:,i) = p.n * p.m * V;
        else
            SY(:,:,i) = p.n * p.m * V;
        end
    end
end

if n_frames > 1
    step = 1/(n_frames-1);
    vx = spline(0:step:1,SX,0:p.dt:1);
    vy = spline(0:step:1,SY,0:p.dt:1);
else
    for i = 1:p.nt
        vx(:,:,i) = SX(:,:,1);
        vy(:,:,i) = SY(:,:,1);
    end
end