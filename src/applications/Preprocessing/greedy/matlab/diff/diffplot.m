function diffplot(I0, I1, wx, wy, fx0, fy0, fx1, fy1, J0, J1)

%{
% Compute \phi^-1 - Id 
[fx0 fy0] = expmap(7, -wx, -wy);

% Compute \phi - Id
[fx1 fy1] = expmap(7, wx, wy);

% Compute I_0 \circ \phi^-1 
J0 = warp(I0, fx0, fy0);

% Compute I_1 \circ \phi 
J1 = warp(I1, fx1, fy1);

%}

% Plot images
subplot(4,2,1);
imagesc(I0); axis image; colormap gray;
title('I_0');

subplot(4,2,2);
imagesc(I1); axis image; colormap gray;
title('I_1');

subplot(4,2,3);
imagesc(J0); axis image; colormap gray;
title('I_0 \circ \phi^{-1}');

subplot(4,2,4);
imagesc(J1); axis image; colormap gray;
title('I_1 \circ \phi');

% Plot grids
subplot(4,2,5);
gridplot(fx0,fy0,5,5);
title('\phi^{-1}');

subplot(4,2,6);
gridplot(fx1,fy1,5,5);
title('\phi');

% Plot W
subplot(4,2,7);
imagesc(wx); axis image; colormap default; caxis([-16 16]);
title('w_x');

subplot(4,2,8);
imagesc(wy); axis image; colormap default; caxis([-16 16]);
title('w_y');
