function [cond] = readCond_3D(cfile,format)
%
% Usage:  [cond] = readCond_3D(cfile,format)
%
%  Reads in Cond3D object, returns as structure
%
%  format = 1 means Mackie's format (default)
%  format = 2 means WS's format

if nargin < 2
    format = 1;
end

if format == 1
    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_mackie3d_model(cfile);
elseif format == 2
    [dx,dy,dz,rho,nzAir,type,origin,rotation] = read_WS3d_model(cfile);
else
    error('Unknown format');
end

grid.dx = dx/1000;
grid.dy = dy/1000;
grid.dz = dz/1000;
grid.Nx = length(dx);
grid.Ny = length(dy);
grid.NzEarth = length(dz);
grid.NzAir = nzAir;
grid.origin = origin/1000;
grid.rotation = rotation;
grid.units = 'km';

if findstr(type,'LOGE')
    cond.paramType = 'LOGE';
    cond.v = - rho;
    cond.AirCond = log(1e-10);
elseif findstr(type,'LOG10')
    cond.paramType = 'LOG10';
    cond.v = - rho;
    cond.AirCond = log10(1e-10);
else
    cond.paramType = 'LINEAR';
    cond.v = 1./rho;
    cond.AirCond = 1e-10;
end
cond.grid = grid;
