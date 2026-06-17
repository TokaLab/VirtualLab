function tokaray = trace_ray_isotropic(tokaray)
% Ray tracing in 2D using the geometrical optics ray equation with scalar
% refractive index
%
%   [Rv, Zv] = trace_ray_isotropic(tokaray)
%
%   Integrates the ray equation (Born & Wolf, Principles of Optics, Cap. 3):
%
%       d/ds ( N * dr/ds ) = grad(N)
%
%   which in normalized form gives:
%
%       dk_hat/ds = (1/N) * ( grad(N) - (grad(N)·k_hat) * k_hat )
%
%   Input:
%       tokaray         struct with fields:
%         .plasma.R     meshgrid R coordinates [m]
%         .plasma.Z     meshgrid Z coordinates [m]
%         .plasma.N0    refractive index on grid [-]
%         .diag.Rin     initial R positions of rays [m]
%         .diag.Zin     initial Z positions of rays [m]
%         .diag.th_in   initial ray angles [rad]
%         .config.ds    integration step [m]
%
%   Output:
%       Rv              ray R positions [n_steps x n_rays]
%       Zv              ray Z positions [n_steps x n_rays]

% unpack
Rg = tokaray.plasma.R;
Zg = tokaray.plasma.Z;
N  = tokaray.plasma.N0;
ds = tokaray.config.ds;

% precompute gradient on grid
[dNdR, dNdZ] = gradient(N, Rg(1,:), Zg(:,1));

% initialize rays
R  = tokaray.diag.Rin;
Z  = tokaray.diag.Zin;
th = tokaray.diag.th_in;

% preallocate
n_steps = tokaray.config.max_steps;
Rv = zeros(n_steps, length(R));
Zv = zeros(n_steps, length(Z));
thv = zeros(n_steps, length(Z));

i = 1;
Rv(i,:) = R;
Zv(i,:) = Z;
thv(i,:) = th;

% initialisation of N for while check
N_i = interp2(Rg, Zg, N, R, Z);

% integration loop
while i < n_steps && ~all(isnan(N_i))

    i = i + 1;

    % interpolate N and gradN at current positions
    N_i    = interp2(Rg, Zg, N,    R, Z);
    dNdR_i = interp2(Rg, Zg, dNdR, R, Z);
    dNdZ_i = interp2(Rg, Zg, dNdZ, R, Z);

    % ray direction unit vector
    k = [cos(th); sin(th)];

    % local gradient
    gradN = [dNdR_i; dNdZ_i];

    % perpendicular component of gradN (drives deflection)
    kdotgradN  = sum(k .* gradN, 1);
    gradN_perp = gradN - kdotgradN .* k;

    % update direction
    k_new = k + (gradN_perp ./ N_i) * ds;
    k_new = k_new ./ vecnorm(k_new);

    % update position
    R  = R + k(1,:) * ds;
    Z  = Z + k(2,:) * ds;
    th = atan2(k_new(2,:), k_new(1,:));

    % record
    Rv(i,:) = R;
    Zv(i,:) = Z;
    thv(i,:) = th;

    % freeze rays that have exited (NaN)
    out = isnan(N_i);
    Rv(i,out)  = Rv(i-1, out);
    Zv(i,out)  = Zv(i-1, out);
    thv(i,out) = thv(i-1,out);
end


% trim to actual length
tokaray.results.R_rt = Rv(1:i-1,:);
tokaray.results.Z_rt = Zv(1:i-1,:);
tokaray.results.th_rt = thv(1:i-1,:);

end