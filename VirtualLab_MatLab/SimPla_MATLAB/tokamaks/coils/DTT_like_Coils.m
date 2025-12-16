function coils = DTT_Coils()
% 
% @Article{en15051702,
% AUTHOR = {Castaldo, Antonio and Albanese, Raffaele and Ambrosino, Roberto and Crisanti, Flavio},
% TITLE = {Plasma Scenarios for the DTT Tokamak with Optimized Poloidal Field Coil Current Waveforms},
% JOURNAL = {Energies},
% VOLUME = {15},
% YEAR = {2022},
% NUMBER = {5},
% ARTICLE-NUMBER = {1702},
% URL = {https://www.mdpi.com/1996-1073/15/5/1702},
% ISSN = {1996-1073},
% ABSTRACT = {In the field of nuclear fusion, the power exhaust problem is still an open issue and represents one of the biggest problems for the realization of a commercial fusion power plant. According to the “European Fusion Roadmap”, a dedicated facility able to investigate possible solutions to heat exhaust is mandatory. For this purpose, the mission of the Divertor Tokamak Test (DTT) tokamak is the study of different solutions for the divertor. This paper presents the plasma scenarios for standard and alternative configurations in DTT. The Single Null scenario is described in detail. The alternative configurations are also presented, showing the good flexibility of the machine.},
% DOI = {10.3390/en15051702}
% }
    
    %% TokaLab Coils
    % Poloidal Field Coils
    PFconfig.names = {"PF1","PF2","PF3","PF4","PF5","PF6"};
    PFconfig.R =  [1.400, 3.080, 4.351,  4.351,  3.080,  1.400];          
    PFconfig.Z =  [2.760, 2.534, 1.015, -1.015, -2.534, -2.760];

    PFconfig.width =  [0.510, 0.279, 0.390, 0.390, 0.279, 0.510];    
    PFconfig.heigth =  [0.590, 0.517, 0.452, 0.452, 0.517, 0.590];
    
    PFconfig.NpixelR = [5, 4, 3, 3, 4, 5];         
    PFconfig.NpixelZ = [6, 5, 4, 4, 5, 6];
    
    
    % Central Solenoid
    CSconfig.names = {"CS3U","CS2U","CS1U","CS1L","CS2L","CS3L"};
    CSconfig.R =  [0.588, 0.588, 0.588, 0.588, 0.588, 0.588];          
    CSconfig.Z =  [2.166, 1.299, 0.433, -0.433, -1.299, -2.166];

    CSconfig.width =  [0.316, 0.316, 0.316, 0.316, 0.316, 0.316];    
    CSconfig.heigth =  [0.788, 0.788, 0.788, 0.788, 0.788, 0.788];
    
    CSconfig.NpixelR = [3, 3, 3, 3, 3, 3];         
    CSconfig.NpixelZ =  [8, 8, 8, 8, 8, 8];

    coils.PFconfig = PFconfig;
    coils.CSconfig = CSconfig;

end