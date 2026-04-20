classdef bremsstrahlung

    properties
        config
    end
     
methods 

    function obj = initialise(obj,method)

        if nargin < 2
            obj.config.method = 1;
        else
            obj.config.method = method;
        end

    end

    function plasma = rad_bremsstrahlung_evalaute(obj,plasma)

        if obj.config.method == 1
            plasma.rad_brem = calculation_m1(plasma);
        else
            % new methods to implement
        end

    end

    %% Bremsstrahlung radiation calculation methods

    function rad_brems = calculation_m1(~,plasma)

        %% Description

        % Computes the bremsstrahlung (free-free) radiation power density in a plasma.
        %
        % This function evaluates a simplified expression for the volumetric
        % bremsstrahlung emission:
        %
        %     P_brems = C * ne * ni * sqrt(Te) * Zeff
        %
        % where:
        %   ne    - electron density [m^-3]
        %   ni    - ion density [m^-3]
        %   Te    - electron temperature (eV)
        %   Zeff  - effective ion charge (dimensionless)
        %   C     - empirical coefficient
        %
        % If Zeff is not provided, it defaults to 1 (pure hydrogen plasma assumption).
        %
        % Input: 
        %   plasma - structure containing the plasma fields   
        %
        % Output:
        %   rad_brems - bremsstrahlung power density [W/m^3]
    
        %% Variables estraction

        ne = plasma.ne; % electron density in [m^{-3}]
        ni = plasma.ni; % ion  density in [m^{-3}]

        Te = plasma.Te;

        if isempty(plasma.Zeff); Zeff = 1; else; Zeff = plasma.Zeff; end

        %% constants
        C = 1.82e-36;

        %% calculation
        rad_brems = C.*ne.*ni.*sqrt(Te).*Zeff;

    end

end



end