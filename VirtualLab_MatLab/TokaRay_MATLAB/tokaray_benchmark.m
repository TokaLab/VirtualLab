%% 

clear; clc;



%%%%%%% benchmark v1 

% 
R = linspace(0,1,100);
Z = linspace(-0.5,0.5,101);
alpha = 0.1;

% plasma definition
[tokaray.plasma.R,tokaray.plasma.Z] = meshgrid(R,Z);
tokaray.plasma.N0 = 1 + alpha*tokaray.plasma.Z; 

% los definition
tokaray.diag.Rin = [0 0 0];
tokaray.diag.Rout = [1 1 1];
tokaray.diag.Zin = [0 0.1 -0.1];
tokaray.diag.Zout = [0 0.3 -0.1];

tokaray.diag.th_in = atan2(tokaray.diag.Zout-tokaray.diag.Zin,...
    tokaray.diag.Rout-tokaray.diag.Rin); 

tokaray.config.ds = 0.001;
tokaray.config.max_steps = 5000;


% code
[Rv, Zv] = trace_ray_isotropic(tokaray);

% analytical solution
R_an = Rv;
Z_an = tokaray.diag.Zin + (sqrt(1+alpha.^2.*R_an.^2)-1)/alpha;

figure(1)
clf
contourf(tokaray.plasma.R,tokaray.plasma.Z,tokaray.plasma.N0,30,'LineStyle','none')
hold on
plot([tokaray.diag.Rin; tokaray.diag.Rout],...
    [tokaray.diag.Zin; tokaray.diag.Zout],'-r')
plot(Rv,Zv,'-b')
plot(R_an,Z_an,'-.k')
grid on
grid minor
xlabel("R [m]")
ylabel("Z [m]")
axis equal



