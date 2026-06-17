function tokaray = import_from_simpla(tokaray,equi)

tokaray.plasma.R = equi.geo.grid.Rg; % [m]
tokaray.plasma.Z = equi.geo.grid.Zg; % [m]

tokaray.plasma.ne = equi.ne; % [m^-3]

tokaray.plasma.Te = equi.Te; % [eV]

tokaray.plasma.Br = equi.Br; % [T]
tokaray.plasma.Bz = equi.Bz; % [T]
tokaray.plasma.Bt = equi.Bt; % [T]

end