function geo = Tokalab_Geometry()
    
    %% Tokalab Geometry
    geo.R0 = 6;
    geo.a = 2;

    geo.grid.kappa_max = 2.2;
    geo.grid.wall_thick = 1.25;

    geo.grid.N_R = 70;
    geo.grid.N_Z = 80;

    geo.wall.R = [3.4, 4.0, 6.0, 8.5, 8.5, 6.0, 4.0, 3.4, 3.4];
    geo.wall.Z = [-4, -5, -5, -3, 3, 5, 5, 4, -4];

end