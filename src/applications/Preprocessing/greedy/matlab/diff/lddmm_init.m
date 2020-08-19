function [vx vy p] = lddmm_init(I0, I1, dt, alpha, gamma, sigma)
% lddmm_init - configure LDDMM problem
% usage:
%   [vx vy p] = lddmm_init(I0, I1, dt, alpha, gamma, sigma)
% parameters:
%   I0, I1      Image 0 and image 1
%   dt          Time step, e.g., 0.05 or 0.1
%   alpha       Weight of the -LBO in metric
%   gamma       Weight of identity in metric
%   sigma       Scaling for image intensity
% returns:
%   vx,vy       Initial time-dependent vector field
%   p           Struct with problem definition, temporaries

% Optional parameters
if nargin < 6, sigma = 1; end;
if nargin < 5, gamma = 1; end;
if nargin < 4, alpha = 10; end;

% Check image dimensions
if size(I0) ~= size(I1)
    error('Images I0 and I1 must have same dimensions!\n');
end
[m n] = size(I0);

% Basic initialization
nt = 1 + 1 / dt;
p = struct(...
    'I0', I0, 'I1', I1,...
    'm',m, 'n', n,...
    'dt',dt,...
    'nt',nt,...
    'alpha',alpha,...
    'gamma',gamma,...
    'sigma',sigma);

% Set up grid
[p.mx p.my] = meshgrid(1:n, 1:m);

% Set FFT of the kernel L LT
p.f_kernel = gamma + 2 * alpha * (2 - cos(2 * pi * (p.mx-1) / m) - cos(2 * pi * (p.my-1) / n));
p.f_kernel_sq = p.f_kernel.^2;

% Image sampling functions
p.image_sample = @(I,x,y)(interp2(I,x,y,'*linear',0));
p.image_sample_gradient = @(I,x,y)(gradient(interp2(I,x,y,'*linear',0)));

% Set up the initial problem
vx = zeros(m,n,nt);
vy = zeros(m,n,nt);