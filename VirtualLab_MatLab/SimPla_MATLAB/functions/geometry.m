classdef geometry
    % Class geometry
    %
    % This class allows to prepare the computational geometry for SimPla
    % (and other module dependent on SimPla, like SynDiag). This class is
    % machine independent since takes the input from the class tokamak
    %
    % Write help geometry.(properties) to see info about each property
    % Write doc geometry to generate the MATLAB documentation for the
    % geometry class

    properties

        R0      % major radius
        a       % minor radius

        R       % major radius (horizontal) coordinate
        Z       % vertical coordinate

        dR      % horizontal step
        dZ      % vertical step

        grid    % structure containing the grid information

        wall    % structure containing the wall information

        operators % structure containing differential operators

    end

    methods

        function obj = import_geometry(obj,tok)
            % Method import_geometry - geometry.import_geometry(tokamak)
            %
            % tokamak is a tokamak class to be initialise and prepared
            % before running this method.
            %
            % The method imports some fundamental information from the
            % tokamak class: R0, a, grid and wall, and store them in the
            % new class geometry

            obj.R0 = tok.R0;
            obj.a = tok.a;
            obj.grid = tok.grid;
            obj.wall = tok.wall;

        end

        function obj = build_geometry(obj)
            % Method build_geometry - geometry.build_geometry()
            %
            % This method builds the computational grid by using R0, a, and
            % the hyperparameters inside grid.
            %
            % See live-scripts or the wiki for details

            R0 = obj.R0;
            a = obj.a;

            wall_thick = obj.grid.wall_thick;
            kappa_max = obj.grid.kappa_max;

            N_R = obj.grid.N_R;
            N_Z = obj.grid.N_Z;

            % horizontal (R) and vertical (Z) coordinates
            R = linspace(R0-a-wall_thick,R0+a+wall_thick,N_R);
            Z = linspace(-kappa_max*a-wall_thick,kappa_max*a+wall_thick,N_Z);

            % generate the grid
            [Rg,Zg] = meshgrid(R,Z);

            % store new information in the class
            obj.grid.Rg = Rg;
            obj.grid.Zg = Zg;

            % prepare differential operator for the grid
            U = utilities;
            [d_dR, d_dZ, d2_dR2, d2_dZ2] = U.differential_operators_fast(obj);

            obj.operators.d_dR = d_dR;
            obj.operators.d_dZ = d_dZ;
            obj.operators.d2_dR2 = d2_dR2;
            obj.operators.d2_dZ2 = d2_dZ2;


            obj.R = R;
            obj.Z = Z;

            obj.dR = R(2)-R(1);
            obj.dZ = Z(2)-Z(1);

        end

        function obj = inside_wall(obj)
            % Method inside_wall - geometry.inside_wall()
            %
            % This method creates a mask of ones for points inside the wall
            % and zeros outside, used several times for both computational
            % and graphics purposes.
            %
            % See live-scripts or the wiki for details

            R = obj.grid.Rg;
            Z = obj.grid.Zg;

            R_wall = obj.wall.R;
            Z_wall = obj.wall.Z;

            inside = inpolygon(R,Z,R_wall,Z_wall);

            obj.wall.inside = inside;
        end

        function [M_wall, indices] = close_to_wall(obj)

            R_wall = obj.wall.R';
            Z_wall = obj.wall.Z';

            R = obj.grid.Rg(:);
            Z = obj.grid.Zg(:);

            % Create a mask for points close to the wall
            d_wall = sqrt((R - R_wall').^2 + (Z - Z_wall').^2);
            [~,indices] = min(d_wall,[],1);

            d = ones(size(R));
            M_wall = spdiags(d', 0, length(d), length(d));
            M_wall = M_wall(indices,:);

        end

        function [M_boundary,indices,ind_bool] = geo_operator(obj)
            % Method inside_wall - geometry.geo_operator()
            %
            % This method generates a mask of ones for the boundary of
            % the geometry and zeros for the rest, used in some cases during
            % Grad-Shavranov equation solver
            %
            % See live-scripts or the wiki for details

            R = obj.grid.Rg;
            Z = obj.grid.Zg;

            % Find boundaries
            ind_bool = (R(:) == R(1,1)) | (Z(:) == Z(1,1)) | (Z(:) == Z(end,1)) |...
                (R(:) == R(1,end)) ;

            indices = find(ind_bool);

            M_boundary = zeros(length(indices),length(ind_bool));

            % Create Boundary Matrix
            for i = 1 : length(indices)
                M_boundary(i,indices(i)) = 1;
            end

            M_boundary = sparse(M_boundary);
        end

        %% Plotting methods

        function plot(obj)
            % Method inside_wall - geometry.plot()
            %
            % Plot the grid points and the wall contours

            plot(obj.grid.Rg(:),obj.grid.Zg(:),'.b')
            hold on
            plot(obj.wall.R, obj.wall.Z, '-k', 'LineWidth', 2) % Plot wall
            axis equal
            grid on
            grid minor
            xlabel('R (Major Radius)');
            ylabel('Z (Vertical Coordinate)');

        end

        function plot_wall(obj)
            % Method inside_wall - geometry.plot_wall()
            %
            % Plot the wall on the open figure

            plot(obj.wall.R, obj.wall.Z, '-k', 'LineWidth', 2) % Plot wall
            grid on
            grid minor
            xlabel('R (Major Radius)');
            ylabel('Z (Vertical Coordinate)');

        end
    end

end