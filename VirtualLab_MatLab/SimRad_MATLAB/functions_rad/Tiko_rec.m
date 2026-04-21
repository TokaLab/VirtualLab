function [Rec] = Tiko_rec(equi,Bolo)

[dpsi_dR, dpsi_dZ] = gradient(equi.psi, equi.geo.dR, equi.geo.dZ);  
gradpsi_mag = sqrt(dpsi_dR.^2 + dpsi_dZ.^2);
nR = dpsi_dR ./ gradpsi_mag;   
nZ = dpsi_dZ ./ gradpsi_mag;   
tR = -nZ; 
tZ =  nR;
% nR_col = nR(:).*equi.geo.wall.inside(:);
% nZ_col = nZ(:).*equi.geo.wall.inside(:);
% tR_col = tR(:).*equi.geo.wall.inside(:);
% tZ_col = tZ(:).*equi.geo.wall.inside(:);

nR_col = nR(:);
nZ_col = nZ(:);
tR_col = tR(:);
tZ_col = tZ(:);

Lx=(equi.geo.operators.d_dR.*tR_col); 
Ly=(equi.geo.operators.d_dZ.*tZ_col); 
Dx=(equi.geo.operators.d_dR.*nR_col); 
Dy=equi.geo.operators.d_dZ.*nZ_col; 




%% Tikhonov

k=[Bolo.Weights';0.0008*(Lx+Ly);  0.0002*(Dx+Dy);equi.geo.wall.inside(:)';];
b=ones(size(k,2),1).*equi.geo.wall.inside(:);
w=zeros(size(k',2),1); w(1:length(Bolo.prj))=Bolo.prj;

for w1=1:5
[b]=lsqr(k,w,1e-12,100,[],[],b);
b(b<0)=0;
end

B = reshape(b,size(equi.geo.grid.Rg));
R0=B; Rec=imgaussfilt(R0,1);
