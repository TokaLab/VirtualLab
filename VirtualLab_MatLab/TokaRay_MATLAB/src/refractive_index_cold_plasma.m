function tokaray = refractive_index_cold_plasma(tokaray)
% REFRACTIVE_INDEX_COLD_PLASMA  Scalar refractive index for unmagnetised plasma
%
%   tokaray = refractive_index_cold_plasma(tokaray)
%
%   Computes the plasma frequency and the scalar refractive index under the
%   cold, unmagnetised plasma approximation:
%
%       omega_pe(x) = sqrt( ne(x) * e^2 / (me * eps0) )
%
%       N^2(x) = 1 - omega_pe^2(x) / omega^2
%
%   This approximation is valid when omega >> omega_ce = eB/me (i.e. the
%   wave frequency is much larger than the electron cyclotron frequency),
%   so the effect of the magnetic field on wave propagation is negligible.
%
%   Regions where N^2 < 0 correspond to cutoff layers where the wave cannot
%   propagate. In this implementation N is set to zero at cutoff (evanescent
%   regions are not modelled).
%
%   Input:
%       tokaray     tokaray object with fields:
%                     .plasma.ne      electron density [m^-3]
%                     .diag.omega     wave angular frequency [rad/s]
%                     .const          physical constants (e_charge, me, eps0)
%
%   Output:
%       tokaray     tokaray object with additional fields:
%                     .plasma.omega_pe    plasma frequency [rad/s]
%                     .plasma.N0          refractive index [-]

ne    = tokaray.plasma.ne;
const = tokaray.const;

% plasma (angular) frequency [rad/s]
tokaray.plasma.omega_pe = sqrt(ne.*const.e_charge.^2./(const.me.*const.eps0));

% scalar refractive index — clamped to zero at cutoff (N^2 < 0)
N2 = 1 - tokaray.plasma.omega_pe.^2./tokaray.diag.omega.^2;
tokaray.plasma.N0 = sqrt(max(N2, 0));

end