classdef geometry

    properties

        R0      % major radius
        a       % minor radius

        R       % major radius (horizontal) coordinate
        Z       % vertical coordinate

        dR      % horizontal step
        dZ      % vertical step

        grid    % structure containing the grid information

        wall    % structure containing the wall information

    end

    methods
        % this function import the geometry from the tokamak class
        function obj = import_geometry(obj,tok)
            obj.R0 = tok.R0;
            obj.a = tok.a;

            obj.grid = tok.grid;

            obj.wall = tok.wall;

        end
        % this function generate the horizontal and vertical coordinates
        % and the grid
        function obj = build_geometry(obj)

            % put object information in local variables
            % (just for improved readability)
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

            obj.R = R;
            obj.Z = Z;

            obj.dR = R(2)-R(1);
            obj.dZ = Z(2)-Z(1);

        end

        % this function calculate a mask for point of the grid inside the
        % wall
        function obj = inside_wall(obj)

            R = obj.grid.Rg;
            Z = obj.grid.Zg;

            R_wall = obj.wall.R;
            Z_wall = obj.wall.Z;

            inside = inpolygon(R,Z,R_wall,Z_wall);

            obj.wall.inside = inside;
        end

        % define the geometry boundary operator
        function [M_boundary,indices,ind_bool] = geo_operator(obj)

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

    end





end