classdef tokaray

% TOKARAY  Ray tracing class for electromagnetic wave propagation in plasma
%
%   Implements geometrical optics ray tracing in 2D and 2.5D for toroidally
%   symmetric plasma. The class acts as a container for plasma data, diagnostic
%   configuration, solver settings, physical constants, and results.
%   Computation is delegated to external functions (one per method) to keep
%   each module self-contained and reusable outside the class.
%
%   Typical workflow from SimPla/SynDiag:
%       tr = tokaray();
%       tr = tr.initialise();
%       tr = tr.from_simpla(equi);              % load plasma from SimPla
%       tr = tr.from_syndiag_ip(IntPol);        % load diagnostic from SynDiag-IP
%       tr = tr.refr_index_cold_plasma();       % compute refractive index
%       tr = tr.raytracing_25D_isotropic();     % trace rays in 2.5D
%       tr = tr.pp_deflection();                % post-process deflection
%
%   For purely poloidal lines of sight (phi_in = phi_out = 0), the 2D solver
%   raytracing_isotropic() is equivalent and slightly faster.
%
%   Properties:
%       plasma      plasma fields on 2D (R,Z) grid:
%                     .R, .Z      spatial coordinates [m]
%                     .ne         electron density [m^-3]
%                     .Te         electron temperature [eV]
%                     .Br/.Bz/.Bt magnetic field components [T]
%                     .omega_pe   plasma frequency [rad/s]
%                     .N0         refractive index [-]
%       diag        diagnostic geometry and wave properties:
%                     .Rin/.Rout/.Zin/.Zout   LOS endpoints [m]
%                     .phi_in/.phi_out         toroidal LOS angles [rad]
%                     .th_in                   poloidal LOS angle [rad]
%                     .k_in                    initial wave vector [3 x n_rays]
%                     .lambda                  wavelength [m]
%                     .omega                   angular frequency [rad/s]
%       config      numerical solver settings:
%                     .ds         integration step [m]
%                     .max_steps  maximum integration steps [-]
%       const       physical constants (see constants.m)
%       results     ray tracing output and post-processing:
%                     .R_rt/.Z_rt/.phi_rt   ray positions [n_steps x n_rays]
%                     .th_rt                poloidal angle [n_steps x n_rays]
%                     .delta_y              linear deflection [m]
%                     .delta_theta          angular deflection [rad]
%
%   References:
%       Born & Wolf, Principles of Optics, Cap. 3
%       Budden, The Propagation of Radio Waves, Cambridge University Press
%       Hutchinson, Principles of Plasma Diagnostics, Cambridge University Press

    properties
        plasma      % plasma fields on 2D (R,Z) grid
        diag        % diagnostic geometry and wave properties
        config      % numerical solver configuration
        const       % physical constants
        results     % ray tracing output and post-processing
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Initialisation function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = initialise(obj)
        % INITIALISE  Load physical constants and set default solver configuration
        %
        %   obj = obj.initialise()
        %
        %   Must be called before any other method.

            % physical constants
            obj.const = constants();

            % default solver configuration
            obj.config.ds        = 0.001;   % integration step [m]
            obj.config.max_steps = 10000;   % maximum integration steps [-]
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Refractive Index Calculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = refr_index_cold_plasma(obj)
        % REFR_INDEX_COLD_PLASMA  Scalar refractive index (unmagnetised limit)
        %
        %   obj = obj.refr_index_cold_plasma()
        %
        %   Valid when omega >> omega_ce. Requires plasma.ne and diag.omega.
        %   Stores omega_pe and N0 in plasma.
        %
        %   See: refractive_index_cold_plasma.m

            obj = refractive_index_cold_plasma(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Ray Tracing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = raytracing_isotropic(obj)
        % RAYTRACING_ISOTROPIC  2D ray tracing in the poloidal plane
        %
        %   obj = obj.raytracing_isotropic()
        %
        %   Integrates the Born & Wolf ray equation in Cartesian-like (R,Z)
        %   coordinates. Valid for lines of sight with no toroidal component
        %   (phi_in = phi_out = 0). Faster than the 2.5D solver.
        %
        %   See: trace_ray_isotropic.m

            obj = trace_ray_isotropic(obj);
        end

        function obj = raytracing_25D_isotropic(obj)
        % RAYTRACING_25D_ISOTROPIC  2.5D ray tracing in cylindrical coordinates
        %
        %   obj = obj.raytracing_25D_isotropic()
        %
        %   Integrates the ray equation in cylindrical (R, phi, Z) coordinates
        %   for toroidally symmetric plasma N = N(R,Z). Handles lines of sight
        %   with a toroidal component (phi_in ~= phi_out). Conserves toroidal
        %   angular momentum L = N*R^2*dphi/ds along each ray.
        %
        %   See: trace_ray_2_5D_isotropic.m

            obj = trace_ray_2_5D_isotropic(obj);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Import functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = from_simpla(obj, equi)
            % FROM_SIMPLA  Import plasma equilibrium from SimPla
            %
            %   obj = obj.from_simpla(equi)
            %
            %   Populates plasma.R, Z, ne, Te, Br, Bz, Bt from a SimPla struct.
            %
            %   See: import_from_simpla.m

            obj = import_from_simpla(obj, equi);
        end

        function obj = from_syndiag_ip(obj, IntPol)
            % FROM_SYNDIAG_IP  Import diagnostic geometry from SynDiag-IP
            %
            %   obj = obj.from_syndiag_ip(IntPol)
            %
            %   Populates diag geometry and wave frequency from a SynDiag-IP struct.
            %   Toroidal angles phi_in and phi_out default to zero (poloidal plane).
            %
            %   See: import_from_syndiag_ip.m

            obj = import_from_syndiag_ip(obj, IntPol);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Post Processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = pp_deflection(obj)
            % PP_DEFLECTION  Post-process ray deflection
            %
            %   obj = obj.pp_deflection()
            %
            %   Computes linear (delta_y) and angular (delta_theta) deflection
            %   of each ray with respect to the unperturbed line of sight.
            %
            %   See: postproc_deflection.m

            obj = postproc_deflection(obj);
        end
        
    end
end