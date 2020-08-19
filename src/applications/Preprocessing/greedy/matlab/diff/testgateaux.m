% Define a random deformation field
rx=padarray(random('norm',0,10,156,156),[50 50]);
ry=padarray(random('norm',0,10,156,156),[50 50]);
wx=20 * imfilter(rx, fspecial('gaussian',80,10));
wy=20 * imfilter(ry, fspecial('gaussian',80,10));
[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));

% Define parameters
[m n] = size(wx);
p.mx = mx;
p.my = my;
p.gamma = 1; p.alpha = 0;
p.A = p.gamma + 2 * p.alpha * (2 - cos(2 * pi * (mx-1) / m) - cos(2 * pi * (my-1) / n));
p.A2 = p.A.^2;
p.sigma = 1;

% Load images
img_c=double(imread('diff_c.png'));
sim_c = imfilter(img_c,fspecial('gaussian',32,4));
img_o=double(imread('diff_o.png'));
sim_o = imfilter(img_o,fspecial('gaussian',32,4));
p.IO = sim_c;
p.I1 = sim_o;

eps = 0.0001;

% Loop over several variations
for i = 1:10
    
    % Define a random variation
    rx=padarray(random('norm',0,10,156,156),[50 50]);
    ry=padarray(random('norm',0,10,156,156),[50 50]);
    hx=20 * imfilter(rx, fspecial('gaussian',80,10));
    hy=20 * imfilter(ry, fspecial('gaussian',80,10));
    
    
    % Compute analytic derivative
    [E, dE] = diff_objective_gateux_deriv(wx, wy, hx, hy, p);
    E2 = diff_objective_gateux_deriv(wx + eps * hx, wy + eps * hy, hx, hy, p);
    E1 = diff_objective_gateux_deriv(wx - eps * hx, wy - eps * hy, hx, hy, p);
    cdE = (E2-E1) / (2*eps);
    
    fprintf('Sol: %12d,  An: %12d   Cd: %12d \n', E, dE, cdE);
end
    