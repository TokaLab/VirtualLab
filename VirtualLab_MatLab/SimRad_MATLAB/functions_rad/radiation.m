classdef radiation

    % radiation description

    properties

        brems

    end

    methods

        %% bremsstrahlung initialisation

        function obj = initialise_brems(obj,method)
            
            if nargin < 2; method = 1; end

            obj.brems = bremsstrahlung();
            obj.brems = obj.brems.initialise(method);
     
        end 

        %% analythical 

       
 
    end
end