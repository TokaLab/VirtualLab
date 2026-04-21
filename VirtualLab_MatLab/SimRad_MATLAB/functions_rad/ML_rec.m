function [R0]=ML_rec(equi,sigma,Bolo)



inside = equi.geo.wall.inside; R=equi.geo.grid.Rg;
Weights=Bolo.Weights';
Weights_t = Weights';
Sn = sum(Weights_t,2) + eps;  
Rv=1000*ones(numel(R),1).*inside(:);
prj=Bolo.prj';

for i = 1:20  
    proj_est = Weights* Rv + eps;
    ratio = prj ./ proj_est;
    backproj = Weights_t * ratio;
    Rv = (Rv ./ Sn) .* backproj;
        
        if sigma > 0 
        R0 = reshape(Rv, size(R));
        R0(inside==0) = 0;
        R0 = imgaussfilt(R0, sigma);
        Rv=R0(:);
        end
   
end



