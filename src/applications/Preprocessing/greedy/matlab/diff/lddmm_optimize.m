function [vxopt vyopt] = lddmm_optimize(vx, vy, p, n_iter)
% lddmm_optimize - run LDDMM optimization
% usage: 
%   [vxopt vyopt] = lddmm_optimize(vx, vy, p)

vxopt = vx;
vyopt = vy;

maxstep = 0.1;

% Line search function
function E = linfun(step)
    E = lddmm_objective_and_gradient(vxopt - step * dedvx,vyopt - step * dedvy, p);
    fprintf('  Line Search : step = %d, val = %d\n', step, E);
end
 
% Iterate
for it = 1:n_iter
    
    % Compute the gradient
    [E dedvx dedvy] = lddmm_objective_and_gradient(vxopt, vyopt, p);
    
    % Print the function value
    fprintf('GradEval %4i: E = %d\n', it, E);   
    
    % Show the results
    lddmm_plot_problem(vxopt, vyopt, p);
    drawnow;
    
    % Search for minimum
    sopt = fminbnd(@linfun, 0, maxstep);
    
    % Break if no improvement
    if sopt == 0, break; end
    
    % Next iteration
    vxopt = vxopt - sopt * dedvx;
    vyopt = vyopt - sopt * dedvy;
end

end