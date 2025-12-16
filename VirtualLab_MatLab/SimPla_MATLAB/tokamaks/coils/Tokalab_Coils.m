function coils = Tokalab_Coils()
    
    %% TokaLab Coils
    
    % Poloidal Field Coils
    PFconfig.names = {"PF1","PF2","PF3","PF4","PF5","PF6"};
    PFconfig.R =  [  4, 8.5, 12, 12, 8.5, 4];          
    PFconfig.Z =  [7.5, 6.5, 3.5, -3.5, -6.5, -7.5];

    PFconfig.width =  [1, 0.65, 0.7, 0.65, 0.8, 1.6];    
    PFconfig.height =  [1, 0.6, 1.1, 1.1, 1, 1];
    
    PFconfig.NpixelR = [4, 3, 3, 3, 3, 6];         
    PFconfig.NpixelZ =  [4, 3, 4, 4, 4, 4];
    
    
    % Central Solenoid
    CSconfig.names = {"CS1","CS2","CS3","CS4","CS5","CS6"};
    CSconfig.R =  [1.7, 1.7, 1.7, 1.7, 1.7, 1.7];          
    CSconfig.Z =  [5, 3, 1, -1, -3, -5];

    CSconfig.width =  [0.75, 0.75, 0.75, 0.75, 0.75, 0.75];    
    CSconfig.heigth =  [2, 2, 2, 2, 2, 2];
    
    CSconfig.NpixelR = [3, 3, 3, 3, 3, 3];         
    CSconfig.NpixelZ =  [8, 8, 8, 8, 8, 8];

    coils.PFconfig = PFconfig;
    coils.CSconfig = CSconfig;

end