function [mvx mvy q] = lddmm_optimize_basic(vx, vy, p, t, n_iter)
% lddmm_optimize - run LDDMM optimization with basic gradient descent
% usage:
%   [vxopt vyopt] = lddmm_optimize(vx, vy, p)
%   use gradient descent to do the optimization

mvx = vx;
mvy = vy;

if n_iter == 0
    [ft0x ft0y ft1x ft1y] = lddmm_integrate_field_semi_lagrangian(vx, vy, p);
    q=struct('ft0x',ft0x,'ft0y',ft0y,'ft1x',ft1x,'ft1y',ft1y);
end;

for ii = 1:n_iter

    % Reparamterize every 10 iterations
    if mod(ii, 10) == 0
        [mvx mvy] = lddmm_reparameterize_vector_field(mvx, mvy, p);
    end
    
    [f mgx mgy q] = lddmm_objective_and_gradient(mvx, mvy, p);
    
    
    mvx = mvx - t * mgx;
    mvy = mvy - t * mgy;

    lddmm_plot_problem_simple(mvx,mvy,p,q);
    getframe;

    if ii > 1
        dcorr = sum(mgx0(:).*mgx(:)+mgy0(:).*mgy(:));
        if  dcorr < 0
            t1 = t;
            t= t * 0.8;

            fprintf(2, 'iter %d/%d, E=%f, dcorr=%f, oscillation: t=%f-->%f \n', ii, n_iter, f,dcorr, t1,t);

            if t < eps
                fprintf(2, 'reach oscillation: t= %f\n', t);
                break;
            end;
        else
            fprintf(2, 'iter %d/%d, E=%f, dcorr=%f, step=%f, \n', ii, n_iter, f,dcorr, t);
        end;
    end;
    

    mgx0 = mgx;
    mgy0 = mgy;

end;
