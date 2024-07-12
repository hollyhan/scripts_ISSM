%function [lambda, mu, BM, YM, PR, VPVS] = lame(vp,vs,rho);
function [lambda, mu] = lame(vp,vs,rho);
% Function to compute elastic constants from seismic velocities and density

% INPUT:
%       vp = P-velocity in km/s
%       vs = S-velocity in km/s
%       rho = density in g/cm^3

% OUTPUT:
%       lambda, mu (Lame's parameters in Pascals)
%       BM (Bulk modulus in Pascals)
%       YM (Young's modulus)
%       PR (Poissonn's ratio)
%       VPVS (the ratio of Vp over Vs) 

% convert velocity & density to MKS units
vp = vp*1000;   % m/s
vs = vs*1000;    % m/s
rho= rho*1000; % kg/m^3

% compute
mu = rho.*vs.*vs;
lambda = vp.^2.*rho - 2*mu;
BM = lambda + 2*mu/3;
YM = (9*BM.*mu) ./ (3*BM + mu);
VPVS = vp./vs;