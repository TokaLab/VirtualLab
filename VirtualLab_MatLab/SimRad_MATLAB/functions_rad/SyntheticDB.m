classdef SyntheticDB
    properties
        N %Number of phantom
        equi
        method
        data
    end
    
    methods
        function obj = SyntheticDB(N, equi,method)
            obj.N = N;
            obj.equi = equi;
            obj.method  = method;
        end
        
        function obj = generate(obj)
                    if obj.method == 'Analytic'
                        rad = radiation();
                        rad = rad.initialise_phantoms(2);
                    for i=1:obj.N

                            sel=randi(4);
                            lambda=rand(1,sel);
                            phantom=zeros(size(obj.equi.geo.grid.Rg));
                                    
                                    for z=1:length(lambda)
                                        rad.analytic.config.I0=1e4+(1e5-1e4)*rand();
                                        rad.analytic.config.mu_p=2*pi*rand();
                                        rad.analytic.config.mu_h=1.005*rand();
                                        rad.analytic.config.std_p=0.1+0.4*rand()*rad.analytic.config.mu_h*pi;
                                        rad.analytic.config.std_h=0.01+0.05*rad.analytic.config.mu_h+0.25*rand()*min(rad.analytic.config.mu_h,0.5);
                                        phantom_1=rad.analytic.calculation_phantom(obj.equi);
                                        phantom=phantom+phantom_1;
                                    end
                                obj.data(i,:)=phantom(:)';   
                    end
                    else
                        return
                    end
        end
       
        
        function save(obj, filename)
            writematrix([obj.data], filename);
        end
  
        function visual_DB(obj)
         
            TP = TokaPlot();
            figura1 = figure(1);
            figura.config.psi_lines = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 1 1.01 1.1];
            figura.config.subplot = [1 1 1];
            figura.config.plot_wall = 1;
            
            figura.config.hold = "on";
            for i=1:size(obj.data,1)
            obj.equi.Rad  = reshape(obj.data(i,:),size(obj.equi.geo.grid.Rg));
            
          
           

            figura1 = TP.PlotField(obj.equi,"Rad",figura1, figura.config);

            end

        end
end

end
