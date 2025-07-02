function geo = TokaPug_Geometry()
    
    %% TokaPug Geometry
    geo.R0 = 5.5;
    geo.a = 2.5;

    geo.grid.kappa_max = 1;
    geo.grid.wall_thick = 1.25;

    geo.grid.N_R = 70;
    geo.grid.N_Z = 80;

    load("TokaPug_wall.mat","TokaPug_geometry")

    geo.wall.R = TokaPug_geometry(:,1)';
    geo.wall.Z = TokaPug_geometry(:,2)';

end