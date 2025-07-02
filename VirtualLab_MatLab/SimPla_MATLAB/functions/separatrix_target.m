classdef separatrix_target

    properties
        R_sep_target
        Z_sep_target

        inside
    end

    methods

        %% build target separatrix
        function obj = build_separatrix(obj,separatrix_config,geo) 
            
            % evaluate separatrix coordinates
            if separatrix_config.method == 1
                obj = build_separatrix_m1(obj,separatrix_config);
            else
                % new method can be added here

            end

            % evaluate grid points inside the separatrix
            obj = inside_separatrix(obj,geo);

        end

        %% Separatrix operators

        function [M_sep,V_sep,ind_sep] = sep_operators(obj,geo)

            R_sep = obj.R_sep_target;
            Z_sep = obj.Z_sep_target;

            R = geo.grid.Rg;
            Z = geo.grid.Zg;

            % we look for all the points closer to the separatrix 
            for i = 1 : length(R_sep)
                d2 = (R-R_sep(i)).^2 + (Z-Z_sep(i)).^2;
                [~,ind_sep(i)] = min(d2,[],'all');
            end

            ind_sep = unique(ind_sep);
            M_sep = zeros(length(ind_sep),size(R(:),1));

            for i = 1 : length(ind_sep)
                M_sep(i,ind_sep(i)) = 1;
                V_sep(i,1) = 0;
            end

        end

        %% Evaluate points inside the separatrix
        function obj = inside_separatrix(obj,geo)
            
            R = geo.grid.Rg;
            Z = geo.grid.Zg;
            
            R_sep = obj.R_sep_target;
            Z_sep = obj.Z_sep_target;

            inside = inpolygon(R,Z,R_sep,Z_sep);

            obj.inside = inside;

        end

        %% this method uses (reference)
        function obj = build_separatrix_m1(obj,separatrix_config)

            % extract parameters
            k1 = separatrix_config.k1;
            k2 = separatrix_config.k2;
            d1 = separatrix_config.d1;
            d2 = separatrix_config.d2;
            gamman1 = separatrix_config.gamma_n_1;
            gamman2 = separatrix_config.gamma_n_2;
            gammap1 = separatrix_config.gamma_p_1;
            gammap2 = separatrix_config.gamma_p_2;

            R0 = separatrix_config.R0;
            Z0 = separatrix_config.Z0;
            a = separatrix_config.a;

            % Inner upper part

            tn1 = (1-d1)/k1*tan(gamman1);

            if tn1 < 0.5
                alpha0n = -(d1-(1+d1)*tn1)/(1-2*tn1);
                alphan = (1-d1)*(1-tn1)/(1-2*tn1);
                betan = k1*(1-tn1)/(sqrt(1-2*tn1));

                thetax = asin(sqrt(1-2*tn1)/(1-tn1));

                theta = linspace(0,thetax,30);

                Rnu = R0 + a*(alpha0n - alphan*cos(theta));
                Znu = a*betan*sin(theta);

            elseif tn1 == 0.5

                zeta = linspace(0,k1,30);
                xi = -1 + (1+d1)/k1^2*zeta^2;

                Rnu = xi*a + R0;
                Znu = zeta*a;

            elseif tn1 > 0.5

                alpha0n = -((1+d1)*tn1-d1)/(2*tn1-1);
                alphan = (1-d1)*(1-tn1)/(2*tn1-1);
                betan = k1*(1-tn1)/(sqrt(2*tn1-1));

                phix = asinh(sqrt(2*tn1-1)/(1-tn1));

                phi = linspace(0,phix,30);

                Rnu = R0 + a*(alpha0n + alphan*cosh(phi));
                Znu = a*betan*sinh(phi);

            elseif tn1 == 1

                zeta = linspace(0,k1,30);
                xi = -1 + (1-d1)/k1 * zeta;

                Rnu = xi*a + R0;
                Znu = zeta*a;

            end

            % Ineer Bottom Part

            tn2 = (1-d2)/k2*tan(gamman2);

            if tn2 < 0.5
                alpha0n = -(d2-(1+d2)*tn2)/(1-2*tn2);
                alphan = (1-d2)*(1-tn2)/(1-2*tn2);
                betan = k2*(1-tn2)/(sqrt(1-2*tn2));

                thetax = asin(sqrt(1-2*tn2)/(1-tn2));

                theta = -thetax:0.01:0;

                Rnl = R0 + a*(alpha0n - alphan*cos(theta));
                Znl = a*betan*sin(theta);

            elseif tn2 == 0.5

                zeta = linspace(-k2,0,30);
                xi = -1 + (1+d2)/k2^2*zeta^2;

                Rnl = xi*a + R0;
                Znl = zeta*a;

            elseif tn2 > 0.5

                alpha0n = -((1+d2)*tn2-d2)/(2*tn2-1);
                alphan = (1-d2)*(1-tn2)/(2*tn2-1);
                betan = k2*(1-tn2)/(sqrt(2*tn2-1));

                phix = asinh(sqrt(2*tn2-1)/(1-tn2));

                phi = linspace(-phix,0,30);

                Rnl = R0 + a*(alpha0n + alphan*cosh(phi));
                Znl = a*betan*sinh(phi);

            elseif tn2 == 1

                zeta = linspace(-k2,0,30);
                xi = -1 + (1-d2)/k2 * zeta;

                Rnl = xi*a + R0;
                Znl = zeta*a;

            end

            % Outer Upper Part

            tp1 = (1+d1)/k1*tan(gammap1);

            if tp1 < 0.5
                alpha0p = -(d1+(1-d1)*tp1)/(1-2*tp1);
                alphap = (1+d1)*(1-tp1)/(1-2*tp1);
                betap = k1*(1-tp1)/(sqrt(1-2*tp1));

                thetax = asin(sqrt(1-2*tp1)/(1-tp1));

                theta = linspace(0,thetax,30);

                Rpu = R0 + a*(alpha0p + alphap*cos(theta));
                Zpu = a*betap*sin(theta);

            elseif tp1 == 0.5

                zeta = linspace(0,k1,30);
                xi = -1 - (1+d1)/k1^2*zeta^2;

                Rpu = xi*a + R0;
                Zpu = zeta*a;

            elseif tp1 > 0.5

                alpha0p = ((1-d1)*tp1+d1)/(2*tp1-1);
                alphap = -(1+d1)*(1-tp1)/(2*tp1-1);
                betap = k1*(1-tp1)/(sqrt(2*tp1-1));

                phix = asinh(sqrt(2*tp1-1)/(1-tp1));

                phi = linspace(0,phix,30);

                Rpu = R0 + a*(alpha0p + alphap*cosh(phi));
                Zpu = a*betap*sinh(phi);

            elseif tn1 == 1

                zeta = linspace(0,k1,30);
                xi = 1 - (1+d1)/k1 * zeta;

                Rpu = xi*a + R0;
                Zpu = zeta*a;

            end

            %Outer Bottom Part

            tp2 = (1+d2)/k2*tan(gammap2);

            if tp2 < 0.5
                alpha0p = -(d2+(1-d2)*tp2)/(1-2*tp2);
                alphap = (1+d2)*(1-tp2)/(1-2*tp2);
                betap = k2*(1-tp2)/(sqrt(1-2*tp2));

                thetax = asin(sqrt(1-2*tp2)/(1-tp2));

                theta = linspace(-thetax,0,30);

                Rpl = R0 + a*(alpha0p + alphap*cos(theta));
                Zpl = a*betap*sin(theta);

            elseif tp2 == 0.5

                zeta = linspace(-k2,0,30);
                xi = -1 - (1+d2)/k2^2*zeta^2;

                Rpl = xi*a + R0;
                Zpl = zeta*a;

            elseif tp2 > 0.5

                alpha0p = ((1-d2)*tp2+d2)/(2*tp2-1);
                alphap = -(1+d2)*(1-tp2)/(2*tp2-1);
                betap = k2*(1-tp2)/(sqrt(2*tp2-1));

                phix = asinh(sqrt(2*tp2-1)/(1-tp2));

                phi = linspace(-phix,0,30);

                Rpl = R0 + a*(alpha0p + alphap*cosh(phi));
                Zpl = a*betap*sinh(phi);

            elseif tp2 == 1

                zeta = linspace(-k2,0,30);
                xi = 1 - (1+d2)/k2 * zeta;

                Rpl = xi*a + R0;
                Zpl = zeta*a;

            end

            % Put Target Separatrix Coordinates togheters
            obj.R_sep_target = [Rnu  flip(Rpu) flip(Rpl) Rnl Rnu(1)];
            obj.Z_sep_target = Z0 + [Znu  flip(Zpu) flip(Zpl) Znl Znu(1)];
           
        end


    end



end