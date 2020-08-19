function [wx wy] = elevolve(img0, img1, gamma, alpha, sigma, epsilon, niter)

% Initialize the first resolution level
maxlev = 1; lfactor = 2^maxlev;
m = size(img0,1) / lfactor;
n = size(img0,2) / lfactor;
wx = zeros(m,n); wy = zeros(m,n);

% Iterate over resolution levels
for lev=1:maxlev

    % Resize the images
    I0 = imresize(imfilter(img0,fspecial('gaussian',256,2^(3-lev))),[m,n]);
    I1 = imresize(imfilter(img1,fspecial('gaussian',256,2^(3-lev))),[m,n]);
    
    % Compute the identity diffeomorphism
    [mx my] = meshgrid(1:n,1:m);

    % Compute the Fourier space multiplier corresponding to Navier operator
    A = gamma + 2 * alpha * (2 - cos(2 * pi * (mx-1) / m) - cos(2 * pi * (my-1) / n));
    A2 = A.*A;

    % Iterative gradient descent
    for it=1:niter

        % Compute \phi^-1 - Id 
        [fx0 fy0] = expmap(7, -wx, -wy);

        % Compute \phi - Id
        [fx1 fy1] = expmap(7, wx, wy);

        % Compute I_0 \circ \phi^-1 
        J0 = interp2(I0, mx + fx0, my + fy0, '*linear', 0);

        % Compute I_1 \circ \phi 
        J1 = interp2(I1, mx + fx1, my + fy1, '*linear', 0);

        % Plotting commands
        diffplot(I0, I1, wx, wy, fx0, fy0, fx1, fy1, J0, J1);
        drawnow;

        % Compute \nabla(I_0 \circ \phi^{-1})
        [Gx0 Gy0] = gradient(J0);

        % Compute \nabla(I_1 \circ \phi)
        [Gx1 Gy1] = gradient(J1);

        % Compute the term on the right hand of L L-adj in (12)
        dIx = (J0 - I1) .* Gx0 - (J1 - I0) .* Gx1;
        dIy = (J0 - I1) .* Gy0 - (J1 - I0) .* Gy1;

        % Apply (LL^{adj})^{-1} to the above
        LLIx = real(ifft2(fft2(dIx) ./ A2));
        LLIy = real(ifft2(fft2(dIy) ./ A2));

        % Compute the update step
        dwx = 2 * (wx - LLIx / sigma^2);
        dwy = 2 * (wy - LLIy / sigma^2);

        % Compute the image match energy
        e_image = (sum(sum((J0 - I1).^2)) + sum(sum((J1 - I0).^2))) / sigma^2;

        % Compute the regularization energy
        e_field = sum(sum(real(ifft2(A .* fft(wx))).^2 + real(ifft2(A .* fft(wy))).^2));

        % Print the energy
        fprintf('Iteration %4i: EField = %12d; EImage = %12d; ETotal = $%12d\n', it, e_field, e_image, e_field+e_image);

        % Update the solution
        wx = wx - epsilon * dwx;
        wy = wy - epsilon * dwy;
    end
    
    % Scale for next resolution level
    m = 2 * m; n = 2 * n;
    wx = 2 * imresize(wx, [m,n],'method','bilinear');
    wy = 2 * imresize(wy, [m,n],'method','bilinear');
    
end