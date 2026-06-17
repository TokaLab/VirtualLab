function tokaray = import_from_syndiag_ip(tokaray,IntPol)

% IMPORT_FROM_SYNDIAG_IP  Import diagnostic geometry from SynDiag
% Inteferometer Polarimeter
%
%   tokaray = import_from_syndiag_ip(tokaray, IntPol)
%
%   Populates the diagnostic fields of a tokaray object from a SynDiag-IP
%   interferometry/polarimetry struct. Imports line-of-sight geometry in
%   cylindrical coordinates (R, phi, Z) and wave frequency.
%
%   The toroidal angles phi_in and phi_out are set to zero by default,
%   corresponding to a line of sight lying in the phi = 0 poloidal plane.
%   For oblique lines of sight, populate diag.phi_in and diag.phi_out
%   manually after calling this function.
%
%   The initial ray angle th_in is the poloidal angle in the (R, Z) plane,
%   computed from the projection of the LOS onto the poloidal plane.
%
%   Input:
%       tokaray     tokaray object (initialised, const must be set)
%       IntPol      SynDiag-IP struct with fields:
%                     .R_in           LOS start R positions [m]
%                     .R_out          LOS end R positions [m]
%                     .Z_in           LOS start Z positions [m]
%                     .Z_out          LOS end Z positions [m]
%                     .config.lambda  wavelength [m]
%
%   Output:
%       tokaray     tokaray object with populated diag fields

% line-of-sight endpoints (transposed to row vectors)
tokaray.diag.Rin = IntPol.R_in';    % [m]
tokaray.diag.Rout = IntPol.R_out';  % [m]
tokaray.diag.Zin = IntPol.Z_in';    % [m]
tokaray.diag.Zout = IntPol.Z_out';  % [m]

% toroidal angles — zero by default (poloidal plane LOS)
tokaray.diag.phi_in = zeros(size(IntPol.Z_out'));   % [rad]
tokaray.diag.phi_out = zeros(size(IntPol.Z_out'));  % [rad]

% poloidal angle of LOS in (R,Z) plane [rad]
tokaray.diag.th_in = atan2(tokaray.diag.Zout-tokaray.diag.Zin,...
    tokaray.diag.Rout-tokaray.diag.Rin); 

% wave frequency from wavelength
tokaray.diag.lambda = IntPol.config.lambda;                     % [m]
tokaray.diag.omega = 2*pi*tokaray.const.c/tokaray.diag.lambda;  % [rad/s]

end