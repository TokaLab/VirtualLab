function tokaray = postproc_deflection(tokaray)

    dR = tokaray.results.R_rt - tokaray.diag.Rin;
    dZ = tokaray.results.Z_rt - tokaray.diag.Zin;
    th = tokaray.diag.th_in;

    tokaray.results.delta_y = dR.*sin(th) - dZ.*cos(th);

    tokaray.results.delta_theta = tokaray.results.th_rt - th;


end