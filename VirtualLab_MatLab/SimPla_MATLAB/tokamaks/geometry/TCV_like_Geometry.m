function geo = TCV_like_Geometry()
    
    %% DTT Geometry
    geo.R0 = 0.88;
    geo.a = 0.25;

    geo.grid.build_method = "wallthick";
    
    geo.grid.R_range = [1.5 4.5];
    geo.grid.Z_range = [-2.4 2.4];

    geo.grid.kappa_max = 2.8;
    geo.grid.wall_thick = 0.2;

    geo.grid.N_R = 70;
    geo.grid.N_Z = 80;

    load("TCV_like_wall.mat","Wall")

    geo.wall.R = Wall(:,1)';
    geo.wall.Z = Wall(:,2)';

end