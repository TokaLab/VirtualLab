classdef equilibrium

    properties

        geo     % structure containing geoemtry information
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

    end

    methods
        
        function obj = import_configuration(obj,geo,config)
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

        function obj = solve_equilibrium(obj,psi)

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

        %% evaluate profiles (MHD and then kinetic)
        function obj = compute_profiles(obj)

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

        %% function to evaluate critical points in standard grid

        function [Opoint,Xpoint] = CriticalPoints(obj,Ip,R,Z,inside_wall,Psi)
            
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

        %% psi post-processing
        function obj = equi_pp(obj)

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

        end

        %% Plotting functions

        % plot target separatrix
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