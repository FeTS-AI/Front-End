function [vxopt vyopt] = lddmm_optimize_lbgfs(vx, vy, p, n_iter)
% lddmm_optimize - run LDDMM optimization
% usage: 
%   [vxopt vyopt] = lddmm_optimize(vx, vy, p)

% Number of variables per vector field
nv = length(vx(:));

% Define internal objective function
function [f G] = my_objective(X, options)
    mvx = reshape(X(1:nv),size(vx));
    mvy = reshape(X(nv+1:end),size(vy));
    if(nargout > 1)
        [f mgx mgy q] = lddmm_objective_and_gradient(mvx, mvy, p);
        G = [mgx(:); mgy(:)];        
        lddmm_plot_problem(mvx,mvy,p,q);
        getframe;
    else
        f = lddmm_objective_and_gradient(mvx, mvy, p);
    end
end

% Perform optimization
my_options = struct('Display','iter','MaxIter',n_iter,'corr',10);
Xopt = minFunc(@my_objective, [vx(:); vy(:)], my_options);
vxopt = reshape(Xopt(1:nv),size(vx));
vyopt = reshape(Xopt(nv+1:end),size(vy));

end

