function [wx wy] = diffopt(I0, I1, gamma, alpha, sigma, epsilon, niter, wx, wy)

% Initialize stationary vector field to zero
[m n] = size(I0);
%wx = zeros(m,n);
%wy = zeros(m,n);

% Compute the identity diffeomorphism
[mx my] = meshgrid(1:n,1:m);

% Compute the Fourier space multiplier corresponding to Navier operator
A = gamma + 2 * alpha * (2 - cos(2 * pi * (mx-1) / m) - cos(2 * pi * (my-1) / n));
A2 = A.*A;

% Construct a parameter vector to pass to objective function
p.I0 = I0;
p.I1 = I1;
p.mx = mx;
p.my = my;
p.A = A;
p.A2 = A2;
p.sigma = sigma;

% Test the objective function
for i = 1:niter

    [f0 dwx dwy] = diff_objective(wx, wy, p);
    
    %Compute derivative in this direction
    [E dE] = diff_objective_gateux_deriv(wx, wy, dwx, dwy, p);
    fprintf('GD in gradient direction: %d\n',dE);
    
    %{
    eps = 0.00001;
    for j = 1:4
        % Create a smooth variation
        rx=padarray(random('norm',0,10,156,156),[50 50]);
        ry=padarray(random('norm',0,10,156,156),[50 50]);
        nux=20 * imfilter(rx, fspecial('gaussian',80,10));
        nuy=20 * imfilter(ry, fspecial('gaussian',80,10));
        
        f1 = diff_objective(wx - nux * eps, wy - nuy * eps, p);
        f2 = diff_objective(wx + nux * eps, wy + nuy * eps, p);
        fprintf('Analytical : %d\n', sum(sum(dwx.*nux)) + sum(sum(dwy.*nuy)));
        fprintf('Frechet    : %d\n', normed_dot_product(dwx, dwy, nux, nuy, p));
        fprintf('CentDiff   : %d\n', (f2 - f1) / (2 * eps));
    end
    %}
    wx = wx - epsilon * dwx;
    wy = wy - epsilon * dwy;
    
    fprintf('Iteration: %3i, Energy: %f\n', i, f0);
    diffplot(p.I0, p.I1, wx, wy); 
    getframe;
    
end

