function [status] = writeCond_2D(cfile,cond,ForPC)
%
% Usage:  [status] = writeCond_2D(cfile,cond)
%
%  writes ModelParam object, provided as structure
%  status is total number of bytes written

nza = cond.grid.Nza;
nz  = cond.grid.Nz;

dy = cond.grid.Dy;
dz = cond.grid.Dz(nza+1:nz);

type = cond.paramType;

if strcmp(type,'LOGE')
    rho = - cond.v;
else
    rho = 1./(cond.v);
end
    
status = write_mackie2d_model(cfile,dy,dz,rho,type);
