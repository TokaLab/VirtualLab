classdef equilibrium
    % Class equilibrium
    %
    % Authors: TokaLab team, 
    % https://github.com/TokaLab/VirtualLab
    % Date: 31/10/2025
    %
    % This class handles most of the computational work in SimPla.
    % It performs the following tasks:
    %   1. Generates the target separatrix based on input from
    %      the tokamak and geometry classes.
    %   2. Solves the plasma equilibrium using a Grad-Shafranov solver.
    %   3. Computes and maps remaining quantities (fields, currents,
    %      pressure, kinetic profiles) based on the selected model.

    properties

        geo     % structure containing geometry information
        config  % structure containing parameters
        separatrix  % separatrix class
        utils   % utils class
        toroidal_curr % methods for evaluate toroidal current
        const   % constant structure
        MHD_prof % class used to evaluate pressure and F2 from psi
        kin_prof % class used to evaluate kinetic profiles (ne, ni, Te, Ti)

        % variables
        psi % poloidal flux [Wb/(2 pi)]
        Psi % poloidal flux [Wb]
        psi_n % normalised poloidal flux

        Jt % toroidal density current
        Jr % radial density current
        Jz % vertical density current

        Bt % toroidal magnetic field
        Br % radial magnetic field
        Bz % vertical magnetic field

        p % pressure
        F2 % F2 field (Grad Shafranov)

        ne % electron density
        ni % ion density

        Te % electron temperature
        Ti % ion temperature

        pe % electron pressure
        pi % ion pressure

        Opoint % O-point
        Xpoint % X-point

        LCFS % Last Close Flux Surface

        psi_ref % Reference psi
        q_psi % profile of q on psi_ref

    end

    methods

        %% initialisation methods

        function obj = import_configuration(obj,geo,config)
            % import_configuration  Import geometry and solver configuration
            %
            %   obj = obj.import_configuration(geo, config)
            %
            %   This method imports the geometry and configuration structures into
            %   the equilibrium object. It also sets default parameters for the
            %   Grad-Shafranov solver, including maximum iterations, tolerances,
            %   update rate, relaxation parameter, and plotting options.
            %
            %   Input:
            %       geo    - Geometry structure (from geometry class)
            %       config - Configuration structure containing solver and tokamak parameters
            %
            %   Output:
            %       obj - Updated equilibrium object with imported geometry and configuration

            obj.geo = geo;
            obj.config = config;

            % equilibrium solver configuration
            obj.config.GSsolver.maxIter = 30;
            obj.config.GSsolver.abs_tol = 0;
            obj.config.GSsolver.rel_tol = 1e-4;
            obj.config.GSsolver.update_rate = 1; % min 0, max 1
            obj.config.GSsolver.Lambda = 0;

            obj.config.GSsolver.Plotting = 1;

        end

        function obj = import_classes(obj)
            % import_classes  Import required classes for equilibrium calculations
            %
            %   obj = obj.import_classes()
            %
            %   This method initializes and imports all auxiliary classes required
            %   for equilibrium computation:
            %       - separatrix_target: defines the target separatrix
            %       - utilities: differential operators and helper functions
            %       - toroidal_current: methods to compute toroidal current
            %       - constants: physical constants
            %       - profile_magnetic: MHD/magnetic profile evaluation
            %       - profile_kinetic: kinetic profile evaluation (ne, ni, Te, Ti)
            %
            %   Output:
            %       obj - Updated equilibrium object with all dependent classes loaded


            % import separatrix class (define the target separatrix)
            obj.separatrix = separatrix_target;

            % import utils class (diff operators)
            obj.utils = utilities;

            % import toroidal current methods
            obj.toroidal_curr = toroidal_current;

            % import constants
            obj.const = constants;

            % import magnetic profiles
            obj.MHD_prof = profile_magnetic;

            % import kinetic profiles
            obj.kin_prof = profile_kinetic;

        end

        %% Grad-Shafranov Solver methods

        function obj = solve_equilibrium(obj,psi)
            % solve_equilibrium  Solve the plasma equilibrium using Grad-Shafranov
            %
            %   obj = obj.solve_equilibrium()
            %   obj = obj.solve_equilibrium(psi)
            %
            %   This method computes the tokamak plasma equilibrium by solving the
            %   Grad-Shafranov equation iteratively. If an initial guess `psi` is not
            %   provided, a first-guess solution is computed automatically.
            %
            %   The method performs the following steps:
            %     1. Builds the target separatrix based on geometry and configuration.
            %     2. Constructs differential operators and boundary/separatrix matrices.
            %     3. Computes an initial guess for the poloidal flux if needed.
            %     4. Iteratively updates the poloidal flux until convergence is reached
            %        according to absolute and relative tolerance criteria.
            %     5. Identifies critical points (O-point, X-point) and normalizes psi.
            %
            %   Input:
            %       psi  - (optional) initial guess of the poloidal flux [Wb/(2pi)]
            %
            %   Output:
            %       obj - Updated equilibrium object with computed fields:
            %               - psi, psi_n (poloidal flux)
            %               - Jt (toroidal current)
            %               - Critical points: Opoint, Xpoint
            %
            %   Notes:
            %       - Solver settings (max iterations, tolerances, update rate, Lambda)
            %         are taken from obj.config.GSsolver.
            %       - Plotting of iterations is enabled if obj.config.GSsolver.Plotting = 1.

            %%

            % psi is the first guess. If it is not given, it is calculated
            % with the method GS_solver.first_guess

            if nargin < 2
                compute_first_guess = 1;
            else
                compute_first_guess = 0;
            end

            %%

            % extract used variable (improved readability)
            R = obj.geo.grid.Rg;
            Z = obj.geo.grid.Zg;
            mu0 = obj.const.mu0;
            Ip = obj.config.toroidal_current.Ip;
            inside_wall = obj.geo.wall.inside;

            % define the separatrix target
            obj.separatrix = obj.separatrix.build_separatrix(obj.config.separatrix,obj.geo);

            % define Grad-Shafranov operator
            [d_dR,~,d2_dR2,d2_dZ2] = obj.utils.differential_operators(obj.geo);
            Delta_star = d2_dR2 - d_dR./R(:) + d2_dZ2;
            Delta_star = sparse(Delta_star);
            clear d_dR d2_dR2 d2_dZ2

            % extract separatrix operator
            [M_sep,V_sep,ind_sep] = obj.separatrix.sep_operators(obj.geo);
            M_sep = sparse(M_sep);

            % extract boundary operator
            [M_boundary, ind_boundary, bool_boundary]= obj.geo.geo_operator();

            % first guess (it applies only if psi is not given as input)
            if compute_first_guess == 1
                % Right-hand term of Grad-Shafranov (normalsied flux equation used)
                Jt = obj.toroidal_curr.Jt_constant(obj.geo,obj.separatrix, obj.config.toroidal_current);
                V_grad = -mu0*R(:).*Jt(:);

                % solve system
                M = [Delta_star; M_sep];
                V = [V_grad; V_sep];

                psi_v = (M'*M)\(M'*V);
                psi = reshape(psi_v,size(R));
            end

            % iterative calculation

            maxIter = obj.config.GSsolver.maxIter;
            abs_tol = obj.config.GSsolver.abs_tol;
            rel_tol = obj.config.GSsolver.rel_tol;
            update_rate = obj.config.GSsolver.update_rate;
            Lambda = obj.config.GSsolver.Lambda;

            convergence = 0;
            iteration = 0;

            % clear figure
            if obj.config.GSsolver.Plotting == 1
                figure(Name="Equi solver")
            end

            while convergence == 0

                iteration = iteration + 1;

                % find critical points (O-point, X-point)
                [Opoint,Xpoint] = obj.CriticalPoints(Ip,R,Z,inside_wall,psi);

                % normalisation of psi
                psi_0 = psi(Opoint);
                psi_b = mean(psi(ind_sep)); % Psi_b = Psi(Xpoint);

                psi_n = (psi-psi_0)./(psi_b-psi_0);

                % update toroidal current given previous psi
                Jt = obj.toroidal_curr.Jt_compute(psi_n,obj.config.toroidal_current,obj.geo,obj.separatrix);

                % New Right-hand term of Grad-Shafranov
                V_grad = -mu0*R(:).*Jt(:);

                % Update Psi on Boundary
                V_boundary = psi(bool_boundary);

                % Compose operator and solution
                M = [Delta_star; M_sep; Lambda*M_boundary];
                V = [V_grad; V_sep; Lambda*V_boundary];

                % solve the system
                psi_v = (M'*M)\(M'*V);
                psi_new = reshape(psi_v,size(R));

                error_abs = mean((psi_new-psi).^2,'all');
                error_rel = error_abs./std(psi,[],'all');

                if error_abs<abs_tol || error_rel<rel_tol || iteration>=maxIter
                    convergence = 1;
                end

                if obj.config.GSsolver.Plotting == 1

                    subplot(1,3,1)
                    hold off
                    contourf(R,Z,psi,30)
                    hold on
                    plot(R(Xpoint),Z(Xpoint),'xr')
                    plot(R(Opoint),Z(Opoint),'or')
                    axis equal
                    xlabel("R [m]")
                    ylabel("z [m]")
                    title("\psi - previous iteration")

                    subplot(1,3,2)
                    hold off
                    contourf(R,Z,psi_new,30)
                    hold on
                    plot(R(Xpoint),Z(Xpoint),'xr')
                    plot(R(Opoint),Z(Opoint),'or')
                    plot(obj.geo.wall.R,obj.geo.wall.Z,'-k','LineWidth',1.2)
                    axis equal
                    xlabel("R [m]")
                    ylabel("z [m]")
                    title("\psi - new iteration")

                    subplot(1,3,3)
                    semilogy(iteration,error_abs,'.b','markersize',16)
                    hold on
                    grid on
                    grid minor
                    xlabel("iteration")
                    ylabel("error [Wb/rad]")

                    drawnow

                end

                % update the psi with update_rate
                psi = update_rate.*psi_new + (1-update_rate)*psi;

            end

            % find critical points (O-point, X-point)
            [Opoint,Xpoint] = obj.CriticalPoints(Ip,R,Z,inside_wall,psi);

            % normalisation of psi
            psi_0 = psi(Opoint);
            psi_b = mean(psi(ind_sep)); % Psi_b = Psi(Xpoint);

            psi_n = (psi-psi_0)./(psi_b-psi_0);

            % save variables
            obj.psi = psi;
            obj.psi_n = psi_n;
            obj.Jt = Jt;

        end

        function obj = solve_equilibrium_dimless(obj,psi)
            %   solve_equilibrium_dimless  Solve the plasma equilibrium using dimensionless Grad-Shafranov
            %
            %   obj = obj.solve_equilibrium_dimless()
            %   obj = obj.solve_equilibrium_dimless(psi)
            %
            %   This method computes the tokamak plasma equilibrium by solving the
            %   Grad-Shafranov equation iteratively. If an initial guess `psi` is not
            %   provided, a first-guess solution is computed automatically.
            %
            %   The method performs the following steps:
            %     1. Builds the target separatrix based on geometry and configuration.
            %     2. Constructs differential operators and boundary/separatrix matrices.
            %     3. Computes an initial guess for the poloidal flux if needed.
            %     4. Iteratively updates the poloidal flux until convergence is reached
            %        according to absolute and relative tolerance criteria.
            %     5. Identifies critical points (O-point, X-point) and normalizes psi.
            %
            %   Input:
            %       psi  - (optional) initial guess of the poloidal flux [Wb/(2pi)]
            %
            %   Output:
            %       obj - Updated equilibrium object with computed fields:
            %               - psi, psi_n (poloidal flux)
            %               - Jt (toroidal current)
            %               - Critical points: Opoint, Xpoint
            %
            %   Notes:
            %       - Solver settings (max iterations, tolerances, update rate, Lambda)
            %         are taken from obj.config.GSsolver.
            %       - Plotting of iterations is enabled if obj.config.GSsolver.Plotting = 1.


            %%

            % psi is the first guess. If it is not given, it is calculated
            % with the method GS_solver.first_guess

            if nargin < 2
                compute_first_guess = 1;
            else
                compute_first_guess = 0;
            end

            %%

            % extract used variable (improved readability)
            R = obj.geo.grid.Rg;
            Z = obj.geo.grid.Zg;
            mu0 = obj.const.mu0;
            Ip = obj.config.toroidal_current.Ip;
            inside_wall = obj.geo.wall.inside;
            R0 = obj.geo.R0;

            % define the separatrix target
            obj.separatrix = obj.separatrix.build_separatrix(obj.config.separatrix,obj.geo);

            % define Grad-Shafranov operator
            [d_dR,~,d2_dR2,d2_dZ2] = obj.utils.differential_operators(obj.geo);
            Delta_star = d2_dR2 - d_dR./R(:) + d2_dZ2;
            Delta_star = sparse(Delta_star);
            clear d_dR d2_dR2 d2_dZ2

            % extract separatrix operator
            [M_sep,V_sep,ind_sep] = obj.separatrix.sep_operators(obj.geo);
            M_sep = sparse(M_sep);

            % extract boundary operator
            [M_boundary, ind_boundary, bool_boundary]= obj.geo.geo_operator();

            % Normalisation factors
            Jc = abs(Ip./R0.^2);
            Psic = abs(mu0.*Jc.*R0.^3);
            DeltaPsic = abs(Psic./R0.^2);
            C_Delta = 1./(DeltaPsic.*length(R(:)));
            C_sep = 1./(Psic.*length(V_sep));

            % first guess (it applies only if psi is not given as input)
            if compute_first_guess == 1
                % Right-hand term of Grad-Shafranov (normalsied flux equation used)
                Jt = obj.toroidal_curr.Jt_constant(obj.geo,obj.separatrix, obj.config.toroidal_current);
                V_grad = -mu0.*R(:).*Jt(:);

                % solve system
                M = [C_Delta.*Delta_star; C_sep.*M_sep];
                V = [C_Delta.*V_grad; C_sep.*V_sep];

                psi_v = (M'*M)\(M'*V);
                psi = reshape(psi_v,size(R));
            end

            % iterative calculation

            maxIter = obj.config.GSsolver.maxIter;
            abs_tol = obj.config.GSsolver.abs_tol;
            rel_tol = obj.config.GSsolver.rel_tol;
            update_rate = obj.config.GSsolver.update_rate;
            Lambda = obj.config.GSsolver.Lambda;

            convergence = 0;
            iteration = 0;

            clear figure
            if obj.config.GSsolver.Plotting == 1
                figure(Name="Equi solver")
            end

            while convergence == 0

                iteration = iteration + 1;

                % find critical points (O-point, X-point)
                [Opoint,Xpoint] = obj.CriticalPoints(Ip,R,Z,inside_wall,psi);

                % normalisation of psi
                psi_0 = psi(Opoint);
                psi_b = mean(psi(ind_sep)); % Psi_b = Psi(Xpoint);

                psi_n = (psi-psi_0)./(psi_b-psi_0);

                % update toroidal current given previous psi
                Jt = obj.toroidal_curr.Jt_compute(psi_n,obj.config.toroidal_current,obj.geo,obj.separatrix);

                % New Right-hand term of Grad-Shafranov
                V_grad = -mu0*R(:).*Jt(:);

                % Update Psi on Boundary
                V_boundary = psi(bool_boundary);

                % Compose operator and solution
                M = [C_Delta.*Delta_star; C_sep.*M_sep; Lambda*M_boundary];
                V = [C_Delta.*V_grad; C_sep.*V_sep; Lambda*V_boundary];

                % solve the system
                psi_v = (M'*M)\(M'*V);
                psi_new = reshape(psi_v,size(R));

                error_abs = mean((psi_new-psi).^2,'all');
                error_rel = error_abs./std(psi,[],'all');

                if error_abs<abs_tol || error_rel<rel_tol || iteration>=maxIter
                    convergence = 1;
                end

                if obj.config.GSsolver.Plotting == 1

                    subplot(1,3,1)
                    hold off
                    contourf(R,Z,psi,30)
                    hold on
                    plot(R(Xpoint),Z(Xpoint),'xr')
                    plot(R(Opoint),Z(Opoint),'or')
                    axis equal
                    xlabel("R [m]")
                    ylabel("z [m]")
                    title("\psi - previous iteration")

                    subplot(1,3,2)
                    hold off
                    contourf(R,Z,psi_new,30)
                    hold on
                    plot(R(Xpoint),Z(Xpoint),'xr')
                    plot(R(Opoint),Z(Opoint),'or')
                    plot(obj.geo.wall.R,obj.geo.wall.Z,'-k','LineWidth',1.2)
                    axis equal
                    xlabel("R [m]")
                    ylabel("z [m]")
                    title("\psi - new iteration")

                    subplot(1,3,3)
                    semilogy(iteration,error_abs,'.b','markersize',16)
                    hold on
                    grid on
                    grid minor
                    xlabel("iteration")
                    ylabel("error [Wb/rad]")

                    drawnow

                end

                % update the psi with update_rate
                psi = update_rate.*psi_new + (1-update_rate)*psi;

            end

            % find critical points (O-point, X-point)
            [Opoint,Xpoint] = obj.CriticalPoints(Ip,R,Z,inside_wall,psi);

            % normalisation of psi
            psi_0 = psi(Opoint);
            psi_b = mean(psi(ind_sep)); % Psi_b = Psi(Xpoint);

            psi_n = (psi-psi_0)./(psi_b-psi_0);

            % save variables
            obj.psi = psi;
            obj.psi_n = psi_n;
            obj.Jt = Jt;

        end

        %% Processing of \psi
        % it finds the last closed surface and critical points O and X points

        % version 1
        function obj = equi_pp(obj)
            % equi_pp  Post-process poloidal flux and compute the Last Closed Flux Surface
            %
            %   obj = obj.equi_pp()
            %
            %   This method post-processes the poloidal flux (psi) after the equilibrium
            %   has been solved. The following steps are performed:
            %     1. Interpolates psi onto a higher-resolution grid for improved accuracy.
            %     2. Determines points inside the wall using geometry information.
            %     3. Calculates critical points (O-point and X-point) on the high-resolution grid.
            %     4. Normalizes psi based on O-point and X-point or separatrix psi values.
            %     5. Identifies the Last Closed Flux Surface (LCFS) by contour analysis.
            %     6. Updates the equilibrium object with normalized flux, LCFS, and critical point coordinates.
            %
            %   Output:
            %       obj - Updated equilibrium object containing:
            %               - psi_n : normalized poloidal flux
            %               - LCFS  : structure with coordinates (R, Z) and inside mask
            %               - Opoint: structure with coordinates of O-point
            %               - Xpoint: structure with coordinates of X-point
            %
            %   Notes:
            %       - Uses target separatrix and wall geometry from obj.separatrix and obj.geo.
            %       - Ensures LCFS and normalized flux are consistent for downstream analysis.

            % variables
            R = obj.geo.grid.Rg;
            Z = obj.geo.grid.Zg;

            R_wall = obj.geo.wall.R;
            Z_wall = obj.geo.wall.Z;

            R_sep_target = obj.separatrix.R_sep_target;
            Z_sep_target = obj.separatrix.Z_sep_target;

            psi = obj.psi;
            Ip = obj.config.toroidal_current.Ip;

            % higher resolution grid
            R_HR = linspace(min(R(:)),max(R(:)),500);
            Z_HR = linspace(min(Z(:)),max(Z(:)),700);
            [R_HR,Z_HR] = meshgrid(R_HR,Z_HR);
            psi_HR = interp2(R,Z,psi,R_HR,Z_HR,"spline");

            % points inside the grid
            inside_wall_HR = inpolygon(R_HR,Z_HR,R_wall,Z_wall);

            % O and X points calculation
            [Opoint,Xpoint] = obj.CriticalPoints(Ip,R_HR, Z_HR,inside_wall_HR,psi_HR);

            Opoint_R = R_HR(Opoint);
            Opoint_Z = Z_HR(Opoint);
            Xpoint_R = R_HR(Xpoint);
            Xpoint_Z = Z_HR(Xpoint);

            %%%%%%%%%%%

            % check if X point is close to target separatrix.
            % if yes, X point is used for normalisation, otherwise we use
            % separatrix psi values
            Xpoint_Sep_distance = min(sqrt((R_sep_target-Xpoint_R).^2 + ...
                (Z_sep_target-Xpoint_Z).^2));

            if Xpoint_Sep_distance < 0.1*obj.geo.a
                psi_O = psi_HR(Opoint);
                psi_X = psi_HR(Xpoint);
            else
                psi_O = psi_HR(Opoint);
                psi_X = mean(interp2(R,Z,psi,R_sep_target,Z_sep_target));
            end

            psi_n = (psi-psi_O)./(psi_X-psi_O);

            %%%%%%%%%%%%%

            % find last closed surface
            level_min = 0.99;
            level_max = 1.01;

            levels = linspace(level_min,level_max,21);

            f = figure;
            lines = contour(R,Z,psi_n,levels);
            close(f);

            k = 0;
            j = 0;
            ind_level = 1;

            while k == 0

                j = j + 1;

                Level = lines(1,ind_level);

                length_level = lines(2,ind_level);

                R_lines = lines(1,ind_level+1:ind_level+length_level);
                Z_lines = lines(2,ind_level+1:ind_level+length_level);

                Closenss = abs(R_lines(1)-R_lines(end)) + abs(Z_lines(1)-Z_lines(end));
                Close = Closenss <= 0.001;

                if Close
                    LCFS.R = R_lines;
                    LCFS.Z = Z_lines;
                end

                ind_level = ind_level + length_level + 1;

                if ind_level >= length(lines(1,:))
                    k = 1;
                end

            end

            LCFS.inside = inpolygon(R,Z,LCFS.R,LCFS.Z);

            obj.LCFS = LCFS;

            obj.Xpoint.R = Xpoint_R;
            obj.Xpoint.Z = Xpoint_Z;
            obj.Opoint.R = Opoint_R;
            obj.Opoint.Z = Opoint_Z;

            obj.psi_n = psi_n;
            obj.Psi = psi*2*pi;

        end

        % version 2
        function obj = equi_pp2(obj)
            % equi_pp2  Post-process poloidal flux using LCFS-based method-
            %           this is a more efficient alternative of equi_pp
            %
            %   obj = obj.equi_pp2()
            %
            %   This method post-processes the poloidal flux after equilibrium
            %   computation by using the Last Closed Flux Surface (LCFS) as reference.
            %   Steps performed:
            %     1. Calls find_LCFS to determine LCFS coordinates.
            %     2. Calculates points inside LCFS at standard resolution.
            %     3. Interpolates psi onto a higher-resolution grid for improved accuracy.
            %     4. Determines critical points (O-point and X-point) based on gradient analysis
            %        and proximity to LCFS.
            %     5. Normalizes psi using O-point and X-point values.
            %     6. Updates the equilibrium object with normalized flux, Psi, and critical points.
            %
            %   Output:
            %       obj - Updated equilibrium object containing:
            %               - psi_n : normalized poloidal flux
            %               - Psi   : poloidal flux [Wb]
            %               - LCFS  : Last Closed Flux Surface (coordinates and inside mask)
            %               - Opoint: coordinates of O-point
            %               - Xpoint: coordinates of X-point
            %
            %   Notes:
            %       - Uses LCFS from obj.LCFS and geometry from obj.geo.
            %       - Ensures consistent normalization of psi for downstream analysis.

            obj = obj.find_LCFS;

            psi = obj.psi;
            R = obj.geo.R;
            Z = obj.geo.Z;
            Ip = obj.config.toroidal_current.Ip;

            % calculate inside LCFS at standard resolution
            obj.LCFS.inside = inpolygon(obj.geo.grid.Rg,obj.geo.grid.Zg,...
                obj.LCFS.R,obj.LCFS.Z);

            % higher resolution grid
            R_HR = linspace(min(R(:)),max(R(:)),500);
            Z_HR = linspace(min(Z(:)),max(Z(:)),700);
            [R_HR,Z_HR] = meshgrid(R_HR,Z_HR);
            Psi = interp2(R,Z,psi,R_HR,Z_HR,"spline");

            inside_LCFS = inpolygon(R_HR,Z_HR,obj.LCFS.R,obj.LCFS.Z);
            inside_wall = inpolygon(R_HR,Z_HR,obj.geo.wall.R,obj.geo.wall.Z);

            % sign correction
            if Ip > 0
                Psi = -Psi;
            end

            % evaluate gradient of poloidal flux
            [dPsidR,dPsidZ] = gradient(Psi);

            R = R_HR;
            Z = Z_HR;

            % find zero values
            gradPsi_2 = (dPsidR./R).^2+(dPsidZ./R).^2;
            ismin = islocalmin(gradPsi_2,1) & islocalmin(gradPsi_2,2) & inside_wall;
            ind = find(ismin);

            % Opoint is defined as the point with the minimum value
            % (negative psi considered) inside the LCFS
            Opoint_ind = find(Psi(ind) == min(Psi(ismin & inside_LCFS)));
            Opoint= ind(Opoint_ind);
            ind(Opoint_ind) = [];

            % Xpoint is defined as closest Xpoint to the LCFS
            [~,Xpoint_ind] = min((obj.LCFS.psi - Psi(ind)).^2);
            Xpoint = ind(Xpoint_ind);

            % if Opoint is not found, geometrical centre is used
            if isempty(Opoint)
                [~,Opoint] = min((R-obj.geo.R0).^2 + (Z-obj.geo.Z0).^2,[],"all");
            end

            % if X point is not found, bottom point of LCSF
            % is used (to be optimised for more generability)
            if isempty(Xpoint)
                [~,Xpoint_bottom] = min(obj.LCFS.Z);
                [~,Xpoint] = min((R-obj.LCFS.R(Xpoint_bottom)).^2 + (Z-obj.LCFS.Z(Xpoint_bottom)).^2,[],"all");
            end

            Opoint_R = R_HR(Opoint);
            Opoint_Z = Z_HR(Opoint);
            Xpoint_R = R_HR(Xpoint);
            Xpoint_Z = Z_HR(Xpoint);

            psi_O = Psi(Opoint);
            psi_X = Psi(Xpoint);

            psi_n = (psi-psi_O)./(psi_X-psi_O);

            obj.Xpoint.R = Xpoint_R;
            obj.Xpoint.Z = Xpoint_Z;
            obj.Opoint.R = Opoint_R;
            obj.Opoint.Z = Opoint_Z;

            obj.psi_n = psi_n;
            obj.Psi = psi*2*pi;

        end

        %% evaluate profiles (MHD and then kinetic)

        function obj = compute_profiles(obj)
            % compute_profiles  Evaluate MHD and kinetic profiles for the equilibrium
            %
            %   obj = obj.compute_profiles()
            %
            %   This method computes the full set of plasma profiles after the equilibrium
            %   has been solved. The following steps are performed:
            %     1. Evaluates pressure (p) and F^2 using the MHD profile class.
            %     2. Computes toroidal magnetic field Bt from F2 and geometry.
            %     3. Computes MHD fields: radial and vertical magnetic fields (Br, Bz)
            %        and corresponding current densities (Jr, Jz).
            %     4. Evaluates kinetic profiles (electron/ion density and temperature)
            %        using the kinetic profile class.
            %     5. Computes electron and ion pressures (pe, pi) from kinetic quantities.
            %
            %   Output:
            %       obj - Updated equilibrium object containing:
            %               - MHD quantities: p, F2, Bt, Br, Bz, Jr, Jz
            %               - Kinetic quantities: ne, ni, Te, Ti, pe, pi
            %
            %   Notes:
            %       - Uses MHD and kinetic profile classes stored in obj.MHD_prof and obj.kin_prof.
            %       - Required for subsequent calculations of q-profile, flux mapping, and diagnostics.

            [p,F2] = obj.MHD_prof.Evaluate_p_F(obj);

            Bt = sign(obj.config.toroidal_current.Bt)*sqrt(F2)./obj.geo.grid.Rg;

            obj.p = p;
            obj.F2 = F2;
            obj.Bt = Bt;

            [Br,Bz,Jr,Jz] = obj.MHD_prof.MHD_fields(obj);

            obj.Br = Br;
            obj.Bz = Bz;
            obj.Jr = Jr;
            obj.Jz = Jz;

            Kinetics = obj.kin_prof.evaluate_profiles(obj);

            obj.ne = Kinetics.ne;
            obj.ni = Kinetics.ni;
            obj.Te = Kinetics.Te;
            obj.Ti = Kinetics.Ti;
            obj.pe = Kinetics.pe;
            obj.pi = Kinetics.pi;

        end

        %% Utilities

        function [Opoint,Xpoint] = CriticalPoints(obj,Ip,R,Z,inside_wall,Psi)
            % CriticalPoints  Identify O-point and X-point in the poloidal flux
            %
            %   [Opoint,Xpoint] = obj.CriticalPoints(Ip,R,Z,inside_wall,Psi)
            %
            %   This method locates the critical points of the poloidal flux:
            %     - O-point: location of the flux minimum (magnetic axis)
            %     - X-point: location of the closest saddle point (magnetic null)
            %
            %   Input:
            %       Ip          - Plasma toroidal current [A]
            %       R, Z        - Grid coordinates (2D arrays)
            %       inside_wall - Logical mask for points inside the wall
            %       Psi         - Poloidal flux on the grid [Wb/(2pi)]
            %
            %   Output:
            %       Opoint - Linear index of the O-point in the grid
            %       Xpoint - Linear index of the X-point in the grid
            %
            %   Notes:
            %       - Psi is inverted for positive Ip to maintain consistent sign.
            %       - Uses local minima of the squared flux gradient to identify points.
            %       - If O-point is not found, the geometric center is used as fallback.
            %       - If X-point is not found, the minimum of the target separatrix is used.

            R_b = obj.separatrix.R_sep_target;
            Z_b = obj.separatrix.Z_sep_target;

            % sign correction
            if Ip > 0
                Psi = -Psi;
            end

            % evaluate gradient of poloidal flux
            [dPsidR,dPsidZ] = gradient(Psi);

            % find zero values
            gradPsi_2 = (dPsidR./R).^2+(dPsidZ./R).^2;
            ismin = islocalmin(gradPsi_2,1) & islocalmin(gradPsi_2,2) & inside_wall;
            ind = find(ismin);

            % Opoint is defined as the point with the minimum value
            % (negative psi considered)
            [~,Opoint_ind] = min(Psi(ismin));
            Opoint = ind(Opoint_ind);
            ind(Opoint_ind) = [];

            % Xpoint is defined as closest Xpoint to the Opoint
            % (alternative methods to be explored (closer to target
            % separatrix?)
            [~,Xpoint_ind] = min((Psi(Opoint) - Psi(ind)).^2);
            Xpoint = ind(Xpoint_ind);

            % if Opoint is not found, geometrical centre is used
            if isempty(Opoint)
                [~,Opoint] = min((R-R0).^2 + (Z).^2,[],"all");
            end

            % if X point is not found, minimum value of target separatrix
            % is used (to be optimised for more generability)
            if isempty(Xpoint)
                [~,Xpoint_boundary] = min(Z_b);
                [~,Xpoint] = min((R-R_b(Xpoint_boundary)).^2 + (Z-Z_b(Xpoint_boundary)).^2,[],"all");
            end

        end

        function obj = find_LCFS(obj)
            % find_LCFS  Identify the Last Closed Flux Surface (LCFS)
            %
            %   obj = obj.find_LCFS()
            %
            %   This method finds the LCFS of the equilibrium using the poloidal flux
            %   distribution (psi) and the machine wall geometry. It performs the
            %   following steps:
            %     1. Computes coarse iso-psi contours inside the wall.
            %     2. Checks for closed contours fully contained inside the wall.
            %     3. Selects the contour with the largest enclosed area as LCFS.
            %     4. Refines the LCFS using a finer set of psi levels around the coarse
            %        estimate.
            %
            %   Output:
            %       obj - Updated equilibrium object with LCFS information:
            %               - LCFS.R   : R coordinates of LCFS
            %               - LCFS.Z   : Z coordinates of LCFS
            %               - LCFS.psi : Poloidal flux value corresponding to LCFS
            %
            %   Notes:
            %       - The method uses contourc to extract iso-psi lines.
            %       - Only closed contours inside the machine wall are considered.
            %       - Refinement ensures a more accurate determination of the LCFS.


            % extract inputs
            psi = obj.psi;
            R = obj.geo.R;
            Z = obj.geo.Z;
            inside = obj.geo.wall.inside;

            % maximum and minimum levels of psi inside the wall
            level_min = min(psi(inside));
            level_max = max(psi(inside));

            % number of levels for first iteration
            n_levels = 30;
            levels_coarse = linspace(level_min,level_max,n_levels);

            % extract iso-psi
            lines = contourc(R,Z,psi,levels_coarse);
            ind_level = 1;
            stop_condition = 0;
            Area = 0;

            while stop_condition == 0

                Level = lines(1,ind_level);
                length_level = lines(2,ind_level);

                R_line = lines(1,ind_level + 1 : ind_level + length_level);
                Z_line = lines(2,ind_level + 1 : ind_level + length_level);

                inside = inpolygon(R_line,Z_line, obj.geo.wall.R ,obj.geo.wall.Z );

                Closness = abs(R_line(1)-R_line(end)) + abs(Z_line(1)-Z_line(end));
                Close = Closness <= 0.001;

                if Close && all(inside)
                    Area_now = abs(polyarea(R_line,Z_line));
                    if Area_now > Area
                        obj.LCFS.R = R_line;
                        obj.LCFS.Z = R_line;
                        obj.LCFS.psi = Level;
                        Area = Area_now;
                    end
                end

                % update new ind level position
                ind_level = ind_level + length_level + 1;

                % stop condition
                if ind_level >= length(lines(1,:))
                    stop_condition = 1;
                end

            end

            % a new iteration is done for psi values close to the
            % previously find LCFS
            i = find(levels_coarse == obj.LCFS.psi);
            n_levels = 30;
            levels_fine = linspace(levels_coarse(i-1),levels_coarse(i+1),n_levels);
            lines = contourc(R,Z,psi,levels_fine);

            ind_level = 1;
            stop_condition = 0;
            Area = 0;

            while stop_condition == 0

                Level = lines(1,ind_level);

                length_level = lines(2,ind_level);

                R_line = lines(1,ind_level + 1 : ind_level + length_level);
                Z_line = lines(2,ind_level + 1 : ind_level + length_level);

                inside = inpolygon(R_line,Z_line, obj.geo.wall.R ,obj.geo.wall.Z );

                Closness = abs(R_line(1)-R_line(end)) + abs(Z_line(1)-Z_line(end));
                Close = Closness <= 0.001;

                if Close && all(inside)
                    Area_now = abs(polyarea(R_line,Z_line));
                    if Area_now > Area
                        obj.LCFS.R = R_line;
                        obj.LCFS.Z = Z_line;
                        obj.LCFS.psi = Level;
                        Area = Area_now;
                    end
                end


                ind_level = ind_level + length_level + 1;

                if ind_level >= length(lines(1,:))
                    stop_condition = 1;
                end

            end

        end

        %% q-profile

        function obj = q_profile(obj)
            % q_profile  Calculate the safety factor (q) profile
            %
            %   obj = obj.q_profile()
            %
            %   This method computes the toroidal safety factor profile q(psi) on the
            %   equilibrium using the poloidal flux (psi), toroidal (Bt) and poloidal
            %   (Bp) magnetic fields. The calculation is performed along iso-psi lines
            %   inside the Last Closed Flux Surface (LCFS).
            %
            %   Steps:
            %     1. Identify points inside the LCFS.
            %     2. Generate a set of psi levels for contour integration.
            %     3. Interpolate Bp and Bt along each iso-psi contour.
            %     4. Compute the line integral to evaluate q for each contour.
            %
            %   Output:
            %       obj.psi_ref - psi levels used for q calculation
            %       obj.q_psi   - safety factor profile corresponding to psi_ref
            %
            %   Notes:
            %       - Points outside the LCFS are assigned psi_n = 2 to exclude them.
            %       - Uses linear interpolation along contour lines to compute the integrals.
            %       - The resulting profile can be used for magnetic stability and confinement analysis.

            % variables for calculation
            psi_n = obj.psi_n;
            Bt = obj.Bt;
            Bp = sqrt(obj.Br.^2 + obj.Bz.^2);

            LCFS = obj.LCFS;

            R = obj.geo.grid.Rg;
            Z = obj.geo.grid.Zg;

            % psi_n = 2 outside from LCFS
            inside_separatrix = inpolygon(R,Z,...
                LCFS.R,LCFS.Z);
            psi_n_c = psi_n;
            psi_n_c(~inside_separatrix) = 2;

            % resolution for integral calculation
            psi_levels = linspace(0.01,1,50);

            % calculation
            for j = 1 : length(psi_levels)

                psi_lines = contour(R,Z,psi_n_c,[-0.1 psi_levels(j)]);

                R_line = psi_lines(1,2:psi_lines(2,1)+1);
                Z_line = psi_lines(2,2:psi_lines(2,1)+1);

                %% integral calculation for q

                Bp_line = interp2(R,Z,Bp,R_line,Z_line);
                Bt_line = interp2(R,Z,abs(Bt),R_line,Z_line);

                R_line_ave = (R_line(1:end-1)+R_line(2:end))/2;
                Bp_line = (Bp_line(1:end-1)+Bp_line(2:end))/2;
                Bt_line = (Bt_line(1:end-1)+Bt_line(2:end))/2;

                dR_line = diff(R_line);
                dZ_line = diff(Z_line);

                dL = sqrt(dR_line.^2+dZ_line.^2);

                q_rec(j) = sum(Bt_line./(R_line_ave.*Bp_line).*dL)/(2*pi);

            end
            close

            obj.psi_ref = psi_levels;
            obj.q_psi = q_rec;

        end

        %% Plotting functions

        function plot_separatrix(obj)

            plot(obj.separatrix.R_sep_target,...
                obj.separatrix.Z_sep_target,'-','linewidth',1.2)
            grid on
            grid minor
            xlabel("R")
            ylabel("Z")

        end

        % plot equilibrium
        function plot_fields(obj,field,equi_lines)

            if nargin < 2
                field = "psi";
                equi_lines = 1;
            elseif nargin < 3
                equi_lines = 1;
            end

            R = obj.geo.grid.Rg;
            Z = obj.geo.grid.Zg;

            F = obj.(field);

            contourf(R,Z,F,30,"LineStyle",'none')
            colorbar()
            hold on
            if equi_lines == 1
                psi_n = obj.psi_n;
                levels = [linspace(0,1,11) 1.01, 1.05 1.1];
                contour(R,Z,psi_n,levels,"r",'linewidth',0.5)
            end
            axis equal
            xlabel("R [m]")
            ylabel("z [m]")
            title(field)
            grid on
            grid minor

        end

    end

end