classdef toroidal_current

    properties

        Jt % toroidal current

    end

    methods

        function Jt = Jt_constant(~,geo,sep,Jt_config)

            dR = geo.dR;
            dZ = geo.dZ;

            inside = sep.inside;
            Ip = Jt_config.Ip;

            Jt_plasma = abs(inside);
            Jt_plasma(~inside) = 0;

            Jt = Jt_plasma.*Ip./sum(Jt_plasma.*dR.*dZ,"all");

        end

        function Jt = Jt_compute(obj,psi_n,Jt_config,geo,sep)

            if Jt_config.method == 1

                Jt = obj.Jt_method_1(psi_n,Jt_config,geo,sep);

            elseif Jt_config.method == 2

                Jt = obj.Jt_method_2(psi_n,Jt_config,geo,sep);

            else
                % here you can implement other methods
            end

        end

        function Jt = Jt_method_1(obj,psi_n,Jt_config,geo,sep)

            % extract variables (improved readability)
            alpha1 = Jt_config.alpha_1;
            alpha2 = Jt_config.alpha_2;
            beta0 = Jt_config.beta_0;
            Ip = Jt_config.Ip;

            dR = geo.dR;
            dZ = geo.dZ;

            R = geo.R;
            R0 = geo.R0;

            % evaluate new Jt given psi (method reference)
            Jt_plasma = (beta0*R/R0 + (1-beta0)*R0./R).*(max(1-psi_n.^alpha1,0)).^alpha2;
            Jt_plasma = Jt_plasma.*sep.inside;
            Jt = Jt_plasma.*Ip./sum(Jt_plasma.*dR.*dZ,"all");

        end

        function Jt = Jt_method_2(obj,psi_n,Jt_config,geo,sep)

            % extract variables (improved readability)
            alpha1 = Jt_config.alpha_1;
            alpha2 = Jt_config.alpha_2;
            beta0 = Jt_config.beta_0;
            Ip = Jt_config.Ip;

            psi_n_peak = Jt_config.psi_n_peak;

            dR = geo.dR;
            dZ = geo.dZ;

            R = geo.R;
            R0 = geo.R0;

            % evaluate new Jt given psi (method reference)
            Jt_plasma = (beta0*R/R0 + (1-beta0)*R0./R).*...
                max(1-((psi_n-psi_n_peak)/(1-psi_n_peak)).^alpha1,0).^alpha2;
            Jt_plasma = Jt_plasma.*sep.inside;
            Jt = Jt_plasma.*Ip./sum(Jt_plasma.*dR.*dZ,"all");

        end



    end

end