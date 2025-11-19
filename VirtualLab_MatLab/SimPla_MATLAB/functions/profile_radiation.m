classdef profile_radiation
   properties


   end


 methods
        function profiles_radiation = evaluate_profiles(obj,equi)

            if equi.config.kinetic.method == 1

                profiles_radiation = obj.profile_radiation_m1(equi);

            else
                % here new methods can be added
            end

        end

        function profiles_radiation = profile_radiation_m1(~,equi)

            
            ne = equi.ne;
            ni = equi.ni;

            Te = equi.Te;
            Ti = equi.Ti;
            
            if isempty(equi.Zeff)
            
                Zeff=1;
            else
            Zeff=equi.Zeff;
            end

            profiles_radiation.Rad = 1.82*1e-36*Zeff*ne.*ni.*sqrt(Te);
            
        end

    end
end

