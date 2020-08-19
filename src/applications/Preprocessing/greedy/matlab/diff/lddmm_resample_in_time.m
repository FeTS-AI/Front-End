function [nvx nvy np] = lddmm_resample_in_time(vx, vy, p, dt_new)
% lddmm_resample_in_time - resample problem and vector field in time
% usage:
%   [nvx nvy np] = lddmm_resample_in_time(vx, vy, p, dt_new)

np = p;
np.dt = dt_new;
np.nt = 1/np.dt + 1;

nvx = shiftdim(interp1(0:p.dt:1, shiftdim(vx,2), 0:np.dt:1),1);
nvy = shiftdim(interp1(0:p.dt:1, shiftdim(vy,2), 0:np.dt:1),1);