
function tokaray = evaluate_direction(tokaray)
% EVALUATE_DIRECTION  Compute initial ray direction in cylindrical coordinates
%
%   tokaray = evaluate_direction(tokaray)
%
%   Computes the unit wave vector k at the ray launch point from the
%   line-of-sight endpoints given in cylindrical coordinates (R, phi, Z).
%
%   The direction is computed by converting endpoints to Cartesian coordinates,
%   taking the difference, normalising, then rotating back to cylindrical:
%
%       kR   =  kx * cos(phi_in) + ky * sin(phi_in)
%       kphi = -kx * sin(phi_in) + ky * cos(phi_in)
%       kZ   =  kz
%
%   Input:
%       tokaray     tokaray object with fields:
%                     .diag.Rin       LOS start R [n_rays] [m]
%                     .diag.Rout      LOS end R   [n_rays] [m]
%                     .diag.Zin       LOS start Z [n_rays] [m]
%                     .diag.Zout      LOS end Z   [n_rays] [m]
%                     .diag.phi_in    LOS start phi [n_rays] [rad]
%                     .diag.phi_out   LOS end phi   [n_rays] [rad]
%
%   Output:
%       tokaray     tokaray object with additional field:
%                     .diag.k_in      initial wave vector [3 x n_rays] [-]
%                                     rows: [kR; kphi; kZ]

% extract LOS endpoints
R_in  = tokaray.diag.Rin;
R_out = tokaray.diag.Rout;
Z_in  = tokaray.diag.Zin;
Z_out = tokaray.diag.Zout;
phi_in  = tokaray.diag.phi_in;
phi_out = tokaray.diag.phi_out;

% convert endpoints from cylindrical to Cartesian
x_in  = R_in  .* cos(phi_in);   y_in  = R_in  .* sin(phi_in);
x_out = R_out .* cos(phi_out);  y_out = R_out .* sin(phi_out);

% direction vector in Cartesian coordinates
dx = x_out - x_in;
dy = y_out - y_in;
dz = Z_out - Z_in;

% normalise to unit vector
nrm = sqrt(dx.^2 + dy.^2 + dz.^2);
dx = dx ./ nrm;
dy = dy ./ nrm;
dz = dz ./ nrm;

% rotate from Cartesian to cylindrical at phi_in
kR_in   =  dx .* cos(phi_in) + dy .* sin(phi_in);
kphi_in = -dx .* sin(phi_in) + dy .* cos(phi_in);
kZ_in   =  dz;

tokaray.diag.k_in = [kR_in; kphi_in; kZ_in];  % [3 x n_rays]
end