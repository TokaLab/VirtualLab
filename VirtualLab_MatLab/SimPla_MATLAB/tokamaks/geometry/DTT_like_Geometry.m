function geo = DTT_like_Geometry()
    
    %% DTT Geometry
    geo.R0 = 2.2;
    geo.a = 0.7;

    geo.grid.kappa_max = 2.1;
    geo.grid.wall_thick = 0.2;

    geo.grid.N_R = 70;
    geo.grid.N_Z = 80;

    load("DTT_like_wall.mat","Wall")

    geo.wall.R = Wall(:,2)';
    geo.wall.Z = Wall(:,1)';

end