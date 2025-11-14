classdef TokaPlot

    properties

    end

    methods

        %% Fields

        function fig = PlotField(~,equi,field,fig,config) 

                disp(nargin)
                
                if nargin < 2
                    disp("missing input")
                    return
                elseif nargin == 2
                    field = "ne";
                    fig = figure()
                    config = struct;
                elseif nargin == 3
                    fig = figure()
                    config = struct;
                elseif nargin == 4
                    config = struct;
                end

                %%
              
                if isfield(config, "subplot")==0
                    config.subplot = [1 1 1];
                end

                figure(fig)

                subplot(config.subplot(1),config.subplot(2), config.subplot(3))
                contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.(field).*equi.geo.wall.inside, 50, "LineStyle", "none")
                hold on

                if isfield(config, "plot_wall")==1 && config.plot_wall == 1

                    x1 = [equi.geo.R(1) equi.geo.R(end) equi.geo.R(end) equi.geo.R(1) equi.geo.R(1)];
                    y1 = [equi.geo.Z(1) equi.geo.Z(1) equi.geo.Z(end) equi.geo.Z(end) equi.geo.Z(1)];
                    patch([x1 equi.geo.wall.R], [y1 equi.geo.wall.Z], [0.75 0.75 0.75], "HandleVisibility","off")
                    plot([2.75 3.4], [-5.65 -4], 'color', [0.75 0.75 0.75], 'LineWidth', 2, "HandleVisibility","off")
                    plot(equi.geo.wall.R, equi.geo.wall.Z, '-k', 'LIneWidth', 1, "HandleVisibility","off")
                    plot([equi.geo.wall.R(1) equi.geo.wall.R(end) equi.geo.wall.R(end) equi.geo.wall.R(1) equi.geo.wall.R(1)], [equi.geo.wall.Z(1) equi.geo.wall.Z(1) equi.geo.wall.Z(end) equi.geo.wall.Z(end) equi.geo.wall.Z(1)], '-k', 'LineWidth', 1, "HandleVisibility","off")

                end
                 if isfield(config, "psi_lines")==1
                    contour(equi.geo.grid.Rg, equi.geo.grid.Zg, equi.psi_n, config.psi_lines, '-w', 'LineWidth', 1.5)
                end

                title(field)
                colorbar()
                colormap("jet")
                axis equal
                xlabel("R [m]")
                ylabel("Z [m]")

        end


        %% Diagnostics

        function fig = PlotDiagnostics(~,diag,fig,config) 

                disp(nargin)
                
                if nargin < 2
                    disp("missing input")
                    return
                elseif nargin == 2
                    field = "ne";
                    fig = figure()
                    config = struct;
                elseif nargin == 3
                    fig = figure()
                    config = struct;
                elseif nargin == 4
                    config = struct;
                end

                %%
              
                if isfield(config, "subplot")==0
                    config.subplot = [1 1 1];
                end

                figure(fig)

                subplot(config.subplot(1),config.subplot(2), config.subplot(3))
                contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.(field).*equi.geo.wall.inside, 50, "LineStyle", "none")
                hold on

                if isfield(config, "plot_wall")==1 && config.plot_wall == 1

                    x1 = [equi.geo.R(1) equi.geo.R(end) equi.geo.R(end) equi.geo.R(1) equi.geo.R(1)];
                    y1 = [equi.geo.Z(1) equi.geo.Z(1) equi.geo.Z(end) equi.geo.Z(end) equi.geo.Z(1)];
                    patch([x1 equi.geo.wall.R], [y1 equi.geo.wall.Z], [0.75 0.75 0.75], "HandleVisibility","off")
                    plot([2.75 3.4], [-5.65 -4], 'color', [0.75 0.75 0.75], 'LineWidth', 2, "HandleVisibility","off")
                    plot(equi.geo.wall.R, equi.geo.wall.Z, '-k', 'LIneWidth', 1, "HandleVisibility","off")
                    plot([equi.geo.wall.R(1) equi.geo.wall.R(end) equi.geo.wall.R(end) equi.geo.wall.R(1) equi.geo.wall.R(1)], [equi.geo.wall.Z(1) equi.geo.wall.Z(1) equi.geo.wall.Z(end) equi.geo.wall.Z(end) equi.geo.wall.Z(1)], '-k', 'LineWidth', 1, "HandleVisibility","off")

                end
                 if isfield(config, "psi_lines")==1
                    contour(equi.geo.grid.Rg, equi.geo.grid.Zg, equi.psi_n, config.psi_lines, '-w', 'LineWidth', 1.5)
                end
                
                end


    end

end