function [fx fy] = expmap_euler(dt, wx, wy)

dx = zeros(size(wx));
dy = zeros(size(wy));
[mx my] = meshgrid(1:size(wx,2),1:size(wx,1));

for t=0:dt:1
    fx = dx + dt * interp2(wx, mx + dx, my + dy, 'linear', 0);
    fy = dy + dt * interp2(wy, mx + dx, my + dy, 'linear', 0);
    dx = fx;
    dy = fy;
    %gridplot(fx,fy,5,5);
    %getframe;
end
    

