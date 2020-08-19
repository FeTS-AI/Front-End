function [wx wy] = elevolve(I0, I1, gamma, alpha, sigma, epsilon, niter)

% Initialize stationary vector field to zero
[m n] = size(I0);
wx = zeros(m,n);
wy = zeros(m,n);

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
    J1 = interp2(I1, mx - fx1, my - fy1, '*linear', 0);
    
    % Plotting commands
    clf;
    
    subplot(3,2,1);
    imshow(I0,'InitialMagnification','fit');
    title('I_0');
    
    subplot(3,2,2);
    imshow(I1,'InitialMagnification','fit');
    title('I_1');
    
    subplot(3,2,3);
    imshow(J0,'InitialMagnification','fit');
    title('I_0 \circ \phi^{-1}');
    
    subplot(3,2,4);
    imshow(J1,'InitialMagnification','fit');
    title('I_1 \circ \phi');
    
    subplot(3,2,5);
    gridplot(fx0,fy0,5,5);
    title('\phi^{-1}');
    
    subplot(3,2,6);
    gridplot(fx1,fy1,5,5);
    title('\phi');
    
    drawnow;

    % Compute \nabla(I_0 \circ \phi^{-1})
    [Gx0 Gy0] = gradient(J0);

    % Compute \nabla(I_1 \circ \phi)
    [Gx1 Gy1] = gradient(J1);

    % Compute the term on the right hand of L L-adj in (12)
    dIx = (J0 - I1) .* Gx0 + (J1 - I0) .* Gx1;
    dIy = (J0 - I1) .* Gy0 + (J1 - I0) .* Gy1;
    
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
