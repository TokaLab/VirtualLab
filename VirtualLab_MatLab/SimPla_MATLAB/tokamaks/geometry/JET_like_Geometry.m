function geo = JET_like_Geometry()
    
    %% JET Geometry
    geo.R0 = 2.85;
    geo.a = 1;

    geo.grid.kappa_max = 1.9;
    
    geo.grid.build_method = "edge";
    geo.grid.R_range = [1.5 4.5];
    geo.grid.Z_range = [-2.4 2.4];

    geo.grid.wall_thick = 0.4;

    geo.grid.N_R = 80;
    geo.grid.N_Z = 90;

    load("JET_like_wall.mat","Wall")

    geo.wall.R = Wall(:,1)';
    geo.wall.Z = Wall(:,2)';

end