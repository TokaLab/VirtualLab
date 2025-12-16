function geo = JET_like_Geometry()
    
    %% JET Geometry
    geo.R0 = 2.85;
    geo.a = 1;

    geo.grid.kappa_max = 1.9;
    geo.grid.wall_thick = 0.2;

    geo.grid.N_R = 70;
    geo.grid.N_Z = 80;

    load("JET_like_wall.mat","Wall")

    geo.wall.R = Wall(:,1)';
    geo.wall.Z = Wall(:,2)';

end