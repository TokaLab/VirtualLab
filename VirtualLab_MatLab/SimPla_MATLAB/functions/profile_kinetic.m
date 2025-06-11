classdef profile_kinetic

    properties


    end

    methods
        function profiles_kinetic = evaluate_profiles(obj,equi)

            if equi.config.kinetic.method == 1

                profiles_kinetic = obj.profile_kinetic_m1(equi);

            else
                % here new methods can be added
            end

        end

        function profiles_kinetic = profile_kinetic_m1(~,equi)

            % variables 
            a1 = equi.config.kinetic.a1;
            a2 = equi.config.kinetic.a2;
            n0 = equi.config.kinetic.n0;
            nsep = equi.config.kinetic.nsep;

            p = equi.p;
            e_charge = equi.const.e_charge;

            psi_n = equi.psi_n;

            % Calculate plasma points 
            inside_LCFS = equi.LCFS.inside;
            inside_wall = equi.geo.wall.inside;

            % Correct Psi
            psi_n_c = psi_n;
            psi_n_c(~inside_LCFS) = 1;
            
            % Evaluate Kinetic Profiles
            n = (n0-nsep).*(1-min(psi_n_c,1).^a1).^a2 + nsep;
            T = p./(2*n.*e_charge);

            n = n.*inside_wall;
            T = T.*inside_wall;

            profiles_kinetic.ne = n;
            profiles_kinetic.ni = n;

            profiles_kinetic.Te = T;
            profiles_kinetic.Ti = T;

            profiles_kinetic.pe = equi.p;
            profiles_kinetic.pi = equi.p;

        end

    end
end

