function A = diffopmat(m, n, alpha, gamma)
% diffopmat - creates sparse matrix for the Navier-Cauchy operator
%   L = \gamma * Id - \alpha \Delta
% usage: A = diffopmat(m, n, alpha, gamma)
%   m,n     Dimensions of the vector field
%   alpha   Weight of the laplacian
%   gamma   Weight of the identity

mn = m * n;

i = reshape(1:mn, [m n]);
i1=circshift(i, [1 0]);
i2=circshift(i, [-1 0]);
i3=circshift(i, [0 1]);
i4=circshift(i, [0 -1]);

r0 = (gamma + 4 * alpha) * ones(mn, 1);
r1 = -alpha * ones(mn,1);

spdat = [...
    i(:) i(:)  r0;
    i(:) i1(:) r1;
    i(:) i2(:) r1;
    i(:) i3(:) r1;
    i(:) i4(:) r1];
    
A=sparse(spdat(:,1), spdat(:,2), spdat(:,3));