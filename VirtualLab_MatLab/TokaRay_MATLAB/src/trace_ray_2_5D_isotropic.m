function tokaray = trace_ray_2_5D_isotropic(tokaray)
% TRACE_RAY_2_5D_ISOTROPIC  Ray tracing in 2.5D with scalar refractive index
%
%   tokaray = trace_ray_2_5D_isotropic(tokaray)
%
%   Integrates the geometrical optics ray equation in cylindrical coordinates
%   (R, phi, Z) for a toroidally symmetric plasma N = N(R, Z).
%
%   The toroidal symmetry of N implies conservation of the toroidal angular
%   momentum along each ray:
%
%       L = N * R^2 * dphi/ds = const
%
%   This reduces the integration to the poloidal plane (R, Z), with an
%   additional centrifugal term in the R equation and phi recovered from L:
%
%       dkR/ds = (1/N) * (dN/dR - kR * (kR*dN/dR + kZ*dN/dZ)) + L^2/(N*R^3)
%       dkZ/ds = (1/N) * (dN/dZ - kZ * (kR*dN/dR + kZ*dN/dZ))
%       dphi/ds = L / (N * R^2)
%
%   kphi is recovered from L at each step; the full 3D unit vector
%   (kR, kphi, kZ) is renormalised after each update.
%
%   Rays that exit the plasma domain (interp2 returns NaN) are frozen at
%   their last valid position. Integration stops when all rays have exited
%   or max_steps is reached.
%
%   Note: when phi_in = phi_out = 0, L = 0 and the centrifugal term vanishes.
%   Results should match trace_ray_isotropic up to floating-point errors
%   (~1e-16) from the normalisation step.
%
%   Input:
%       tokaray     tokaray object with fields:
%                     .plasma.R           R meshgrid [m]
%                     .plasma.Z           Z meshgrid [m]
%                     .plasma.N0          refractive index on grid [-]
%                     .diag.Rin           initial R positions [m]
%                     .diag.Zin           initial Z positions [m]
%                     .diag.phi_in        initial phi positions [rad]
%                     .diag.th_in         initial poloidal angle [rad]
%                     .diag.Rout          LOS end R (for evaluate_direction) [m]
%                     .diag.Zout          LOS end Z (for evaluate_direction) [m]
%                     .diag.phi_out       LOS end phi (for evaluate_direction) [rad]
%                     .config.ds          integration step [m]
%                     .config.max_steps   maximum number of steps [-]
%
%   Output:
%       tokaray     tokaray object with additional fields:
%                     .results.R_rt       R positions   [n_steps x n_rays] [m]
%                     .results.Z_rt       Z positions   [n_steps x n_rays] [m]
%                     .results.phi_rt     phi positions [n_steps x n_rays] [rad]
%                     .results.th_rt      poloidal angle [n_steps x n_rays] [rad]
%
%   References:
%       Born & Wolf, Principles of Optics, Cap. 3

% unpack plasma grid and refractive index
Rg = tokaray.plasma.R;
Zg = tokaray.plasma.Z;
N  = tokaray.plasma.N0;
ds = tokaray.config.ds;

% precompute gradient of N on full grid
[dNdR, dNdZ] = gradient(N, Rg(1,:), Zg(:,1));

% initialise ray positions
R  = tokaray.diag.Rin;
Z  = tokaray.diag.Zin;
phi = tokaray.diag.phi_in;
th = tokaray.diag.th_in;

% compute initial wave vector components from LOS geometry
tokaray = evaluate_direction(tokaray);
k_in = tokaray.diag.k_in;
kR = k_in(1,:);
kphi = k_in(2,:);
kZ = k_in(3,:);

% preallocate output arrays [max_steps x n_rays]
n_steps = tokaray.config.max_steps;
n_rays  = length(R);
Rv   = zeros(n_steps, n_rays);
Zv   = zeros(n_steps, n_rays);
phiv = zeros(n_steps, n_rays);
thv  = zeros(n_steps, n_rays);

% compute conserved toroidal angular momentum: L = N * R^2 * kphi
N_0 = interp2(Rg, Zg, N, R, Z);
L   = N_0 .* R.^2 .* kphi;     % [1 x n_rays] — constant along each ray

% record initial state
i = 1;
Rv(i,:)   = R;
Zv(i,:)   = Z;
phiv(i,:) = phi;
thv(i,:)  = th;

% initial interpolation needed for while condition
N_i = interp2(Rg, Zg, N, R, Z);

% integration loop — stops when all rays exit domain or max_steps is reached
while i < n_steps && ~all(isnan(N_i))

    i = i + 1;

    % interpolate N and its gradient at current ray positions
    N_i    = interp2(Rg, Zg, N,    R, Z);
    dNdR_i = interp2(Rg, Zg, dNdR, R, Z);
    dNdZ_i = interp2(Rg, Zg, dNdZ, R, Z);

    % projection of grad(N) along poloidal ray direction
    kdotgradN = kR.*dNdR_i + kZ.*dNdZ_i;

    % poloidal ray equations in cylindrical coordinates
    % centrifugal term L^2/(N*R^3) accounts for toroidal momentum
    dkdR = (dNdR_i - kR.*kdotgradN) ./ N_i + L.^2 ./ (N_i .* R.^3);
    dkdZ = (dNdZ_i - kZ.*kdotgradN) ./ N_i;

    % update poloidal components (Euler step)
    kR = kR + dkdR * ds;
    kZ = kZ + dkdZ * ds;

    % recover kphi from conservation of L, then renormalise full 3D vector
    kphi  = L ./ (N_i .* R.^2);
    knorm = sqrt(kR.^2 + (R.*kphi).^2 + kZ.^2);
    kR    = kR   ./ knorm;
    kZ    = kZ   ./ knorm;
    kphi  = kphi ./ knorm;

    % update positions
    R   = R   + kR  .* ds;
    Z   = Z   + kZ  .* ds;
    phi = phi + L ./ (N_i .* R.^2) .* ds;  % dphi/ds = L/(N*R^2)
    th  = atan2(kZ, kR);                    % poloidal angle in (R,Z) plane


    th = atan2(kZ, kR);

    % record
    Rv(i,:)   = R;
    Zv(i,:)   = Z;
    phiv(i,:) = phi;
    thv(i,:)  = th;

    % freeze rays that have exited (NaN)
    out = isnan(N_i);
    Rv(i,out)  = Rv(i-1, out);
    Zv(i,out)  = Zv(i-1, out);
    thv(i,out) = thv(i-1,out);
    phiv(i,out) = phiv(i-1,out);
end


% trim preallocated arrays to actual integration length
tokaray.results.R_rt = Rv(1:i-1,:);
tokaray.results.Z_rt = Zv(1:i-1,:);
tokaray.results.th_rt = thv(1:i-1,:);
tokaray.results.phi_rt = phiv(1:i-1,:);

end