clear; close all;
%%
load("DB1.mat")
box.R= [2.75 2.75 9.25 9.25]; box.Z= [-5.65 5.65 5.65 -5.65];
equi=DB1.equi;

Bolo = Diag_Bolo();
Bolo  = Bolo.Upload(1);

for i=1:(DB1.N)
equi.Rad = reshape(DB1.data(i,:), size(equi.geo.grid.Rg));
Bolo = Bolo.measure(equi);
R0 = ML_rec(equi,1,Bolo);
R1 = Tiko_rec(equi,Bolo);
Phantom_plot= equi.Rad;
Phantom_plot(~equi.geo.wall.inside)=NaN;
Rec_tiko_plot = R1; Rec_tiko_plot(~equi.geo.wall.inside)=NaN;
Rec_ml_plot = R0; Rec_ml_plot(~equi.geo.wall.inside)=NaN;

f=figure(1)
f.Position= [152	271	1384	550];
tiledlayout('horizontal')

subplot(1,3,1)
hold off
fill(box.R,box.Z,[0.75 0.75 0.75]);hold on
fill(equi.geo.wall.R,equi.geo.wall.Z,[0.2422 0.1504 0.6603]);  hold on
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,Phantom_plot,'LineStyle','none');
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,40,'k','LineWidth',0.5)
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',2)
plot(equi.geo.wall.R,equi.geo.wall.Z,'k','LineWidth',4); 
%plot(equi.LCFS.R,equi.LCFS.Z,'r--','LineWidth',4)
axis equal
C2=colorbar()

ylabel(C2,'[a.u]')
xlabel('R [m]')
ylabel('Z [m]')
title('Phantom')



subplot(1,3,2)
hold off
fill(box.R,box.Z,[0.75 0.75 0.75]);hold on
fill(equi.geo.wall.R,equi.geo.wall.Z,[0.2422 0.1504 0.6603]);  hold on
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,Rec_tiko_plot,'LineStyle','none');
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,20,'k','LineWidth',0.5)
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',2)
plot(equi.geo.wall.R,equi.geo.wall.Z,'k','LineWidth',4); 
%plot(equi.LCFS.R,equi.LCFS.Z,'r--','LineWidth',4)
axis equal
C2=colorbar()
%caxis([0 1e5])
ylabel(C2,'W/m^3')
xlabel('R [m]')
ylabel('Z [m]')
title('Reconstruction Tiko')

subplot(1,3,3)
hold off
fill(box.R,box.Z,[0.75 0.75 0.75]);hold on
fill(equi.geo.wall.R,equi.geo.wall.Z,[0.2422 0.1504 0.6603]);  hold on
contourf(equi.geo.grid.Rg,equi.geo.grid.Zg,Rec_ml_plot,'LineStyle','none');
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,20,'k','LineWidth',0.5)
contour(equi.geo.grid.Rg,equi.geo.grid.Zg,equi.psi_n,[-0.1 1],...
    'r','LineWidth',2)
plot(equi.geo.wall.R,equi.geo.wall.Z,'k','LineWidth',4); 
%plot(equi.LCFS.R,equi.LCFS.Z,'r--','LineWidth',4)
axis equal
C2=colorbar()
%caxis([0 1e5])
ylabel(C2,'W/m^3')
xlabel('R [m]')
ylabel('Z [m]')
title('Reconstruction ML')
% 
% 
sgtitle(num2str(i))
drawnow








end
