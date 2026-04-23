clear; clc;

%% Upload the dataset configuration

db = Benchmarking_Tomography_v1();

%% generation

db_gen = db_gen_radiation;
db = db_gen.generate(db);

%% save 

db_gen.save_db(db,"riccardo")

%% plot

figure(1)
clf
subplot(1,3,1)
g = 1; k = 1;
contourf(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    reshape(db{g}.data_rad(k,:),size(db{g}.equi.psi_n)),30,"LineStyle",'none')
hold on
contour(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    db{g}.equi.psi_n,[0 0.25 0.5 0.75 0.9 0.95 0.99 1 1.01 1.05],...
    '-w','LineWidth',0.5)
plot(db{g}.equi.geo.wall.R,db{g}.equi.geo.wall.Z,'-k','LineWidth',2)
axis equal
xlabel("R [m]")
ylabel("Z [m]")
grid on
grid minor
colorbar()

subplot(1,3,2)
g = 2; k = 3;
contourf(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    reshape(db{g}.data_rad(k,:),size(db{g}.equi.psi_n)),30,"LineStyle",'none')
hold on
contour(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    db{g}.equi.psi_n,[0 0.25 0.5 0.75 0.9 0.95 0.99 1 1.01 1.05],...
    '-w','LineWidth',0.5)
plot(db{g}.equi.geo.wall.R,db{g}.equi.geo.wall.Z,'-k','LineWidth',2)
axis equal
xlabel("R [m]")
ylabel("Z [m]")
grid on
grid minor
colorbar()

subplot(1,3,3)
g = 3; k = 3;
contourf(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    reshape(db{g}.data_rad(k,:),size(db{g}.equi.psi_n)),30,"LineStyle",'none')
hold on
contour(db{g}.equi.geo.grid.Rg,db{g}.equi.geo.grid.Zg,...
    db{g}.equi.psi_n,[0 0.25 0.5 0.75 0.9 0.95 0.99 1 1.01 1.05],...
    '-w','LineWidth',0.5)
plot(db{g}.equi.geo.wall.R,db{g}.equi.geo.wall.Z,'-k','LineWidth',2)
axis equal
xlabel("R [m]")
ylabel("Z [m]")
grid on
grid minor
colorbar()


%%%%%% Cose che devo fare

% equilibrio configurabile

%%%%%% Domande per Wyss

% generazione randomica, numeri statici? 


