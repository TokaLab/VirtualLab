function coils = DTT_like_Coils()
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
    

% @article{AMBROSINO2023113714,
% title = {Conceptual design of the DTT in-vessel equatorial coils},
% journal = {Fusion Engineering and Design},
% volume = {194},
% pages = {113714},
% year = {2023},
% issn = {0920-3796},
% doi = {https://doi.org/10.1016/j.fusengdes.2023.113714},
% url = {https://www.sciencedirect.com/science/article/pii/S0920379623002971},
% author = {R. Ambrosino and R. Albanese and E. Acampora and A. Castaldo and F. Crisanti and R. Iervolino and A. Lampasi},
% keywords = {DTT tokamak, Vertical stability control, IN-vessel coil, Magnetic control system},
% abstract = {For the next generation fusion devices, the presence of in-vessel coils is fundamental to deal with elongated plasmas and alternative configurations, highly performant, yet critically sensitive to the vertical unstable mode. A pair of equatorial in-vessel coils connected in anti-series is sufficient to deal with the vertical control problem. In case of independent equatorial coils, they can also be used to produce a fast radial control action able to preserve plasma facing components during fast plasma transients. The main criticality of independent in-vessel coils is related to the induced current during disruptive events, making the coils extremely vulnerable due to overcurrents and mechanical stress. In this paper, we present the solution proposed for the conceptual design of the in-vessel equatorial coils for the Divertor Tokamak Test fusion device. This solution is able to guarantee radial and vertical control performance and robustness in case of disruption. Simulations are provided to show the efficiency of the design in case of disruptions and the effectiveness of the closed-loop control actions.}
% }


    %% TokaLab Coils
    % Poloidal Field Coils
    PFconfig.names = {"PF1","PF2","PF3","PF4","PF5","PF6","IVC1","IVC2","IVC3","IVC4","VSL","VSU"};
    PFconfig.R =  [1.400, 3.080, 4.351,  4.351,  3.080,  1.400, 1.45, 1.70, 2.25, 2.85, 3.108, 3.108];          
    PFconfig.Z =  [2.760, 2.534, 1.015, -1.015, -2.534, -2.760, -1.15,-1.75,-1.65,-1.00, -0.553, 0.553];

    PFconfig.width =  [0.510, 0.279, 0.390, 0.390, 0.279, 0.510, 0.07, 0.07, 0.07, 0.07, 0.104, 0.104];    
    PFconfig.heigth = [0.590, 0.517, 0.452, 0.452, 0.517, 0.590, 0.07, 0.07, 0.07, 0.07, 0.130, 0.130];
    
    PFconfig.NpixelR = [5, 3, 4, 4, 3, 5, 1, 1, 1, 1, 1, 1];         
    PFconfig.NpixelZ = [6, 5, 4, 4, 5, 6, 1, 1, 1, 1, 2, 2];
    
    
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