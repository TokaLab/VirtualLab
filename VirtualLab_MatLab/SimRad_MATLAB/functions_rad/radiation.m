classdef radiation

    % radiation description

    properties
        brems
        analytic
    end

    methods

        %% bremsstrahlung initialisation

        function obj = initialise_brems(obj,method)
            
            if nargin < 2; method = 1; end

            obj.brems = bremsstrahlung();
            obj.brems = obj.brems.initialise(method);
     
        end 

        %% analythical 
        
        function obj = initialise_phantoms(obj,method)
            
            if nargin < 2; method = 1; end

            obj.analytic = phantom();
            obj.analytic = obj.analytic.initialise(method);
     
        end 

       
 
    end
end