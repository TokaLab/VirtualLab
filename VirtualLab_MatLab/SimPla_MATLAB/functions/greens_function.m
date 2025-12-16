classdef greens_function

    properties

        Green_psi
        Green_Br
        Green_Bz


    end

    methods

        %% Plasma
        function [Green_psi, Green_Br, Green_Bz, ind_J] = GreensFunctionPlasma_method1(~,equi,ind)

            % extract important variables
            mu0 = equi.const.mu0;

            ind_J = find(abs(equi.Jt)>0);

            R = equi.geo.grid.Rg(ind_J);
            Z = equi.geo.grid.Zg(ind_J);

            Rc = equi.geo.grid.Rg(ind)';
            Zc =  equi.geo.grid.Zg(ind)';


            % Safety cutoff values
            eps_val = 1e-9;       % small offset for zero division
            k_max = 0.9999;       % cutoff for elliptic integral singularity

            % calculates Greens function
            a = sqrt((R+Rc).^2+(Z-Zc).^2);
            a(a < eps_val) = eps_val; % avoid zero distance


            k2 = (4 .* R .* Rc) ./ (a.^2);
            k2(k2 > 1) = 1; % theoretical limit
            k = sqrt(k2);
            k(k > k_max) = k_max; % avoid singularity near k=1
            [K, E] = ellipke(k);


            % Avoid division by zero in T2
            denom = a.^3 .* (1 - k);
            denom(denom < eps_val) = eps_val;

            T1 = K ./ a;
            T2 = E ./ denom;


            Green_psi = (mu0./(2*pi)).*a.*((1-k/2).*K - E);
            Green_Br = (mu0/(2*pi)).*((Z-Zc).*(T2.*(Z-Zc).^2 + R.^2 + Rc.^2) - T1)./(R);
            Green_Bz = (mu0/(2*pi)).*(T1 + (Rc.^2 - R.^2 - (Z-Zc).^2).*T2);

        end

        %% Coils

        function [Green_psi, Green_Br, Green_Bz] = GreensFunctionCoils_method1(~,equi,coils)

            % extract important variables
            mu0 = equi.const.mu0;

            Rg = equi.geo.grid.Rg(:);
            Zg = equi.geo.grid.Zg(:);
            nPoints = numel(Rg);

            % Safety cutoff values
            eps_val = 1e-9;       % small offset for zero division
            k_max = 0.9999;       % cutoff for elliptic integral singularity


            coilNames = fieldnames(coils.system);
            nCoils = numel(coilNames);

            % Initialise output matrices
            Green_psi = zeros(nPoints, nCoils);
            Green_Br  = zeros(nPoints, nCoils);
            Green_Bz  = zeros(nPoints, nCoils);

            % Loop on each coil
            for c = 1:nCoils
                coil = coils.system.(coilNames{c});

                [Rc_mesh, Zc_mesh] = meshgrid(coil.R, coil.Z);
                Rc_pixels = Rc_mesh(:);
                Zc_pixels = Zc_mesh(:);

                nPixels = numel(Rc_pixels);

                psi_sum = zeros(nPoints,1);
                Br_sum  = zeros(nPoints,1);
                Bz_sum  = zeros(nPoints,1);

                % Loop on pixels
                for p = 1:nPixels
                    a = sqrt((Rg + Rc_pixels(p)).^2 + (Zg - Zc_pixels(p)).^2);
                    a(a < eps_val) = eps_val; % avoid zero distance

                    k2 = (4 .* Rg .* Rc_pixels(p)) ./ (a.^2);
                    k2(k2 > 1) = 1; % theoretical limit
                    k = sqrt(k2);
                    k = sqrt(k2);
                    k(k > k_max) = k_max; % avoid singularity near k=1

                    [K, E] = ellipke(k);

                    % Avoid division by zero in T2
                    denom = (a.^3 .* (1 - sqrt(k)));
                    denom(denom < eps_val) = eps_val;

                    T1 = K ./ a;
                    T2 = E ./ denom;

                    psi_pixel = (mu0/(2*pi)) .* a .* ((1 - sqrt(k)/2).*K - E);
                    Br_pixel  = (mu0/(2*pi)) .* ((Zg - Zc_pixels(p)) .* ...
                        (T2 .* ((Zg - Zc_pixels(p)).^2 + Rg.^2 + Rc_pixels(p).^2) - T1)) ./ max(Rg, eps_val);
                    Bz_pixel  = (mu0/(2*pi)) .* (T1 + (Rc_pixels(p).^2 - Rg.^2 - (Zg - Zc_pixels(p)).^2) .* T2);

                    psi_sum = psi_sum + psi_pixel;
                    Br_sum  = Br_sum  + Br_pixel;
                    Bz_sum  = Bz_sum  + Bz_pixel;
                end

                Green_psi(:,c) = psi_sum./nPixels;
                Green_Br(:,c)  = Br_sum./nPixels;
                Green_Bz(:,c)  = Bz_sum./nPixels;
            end
        end



    end
end

