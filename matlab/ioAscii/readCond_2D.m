function [cond] = readCond_2D(cfile,PLOT)
%
% Usage:  [cond] = readCond_2D(cfile,PLOT)
%
%  Reads in Cond2D object, returns as structure
%  Second argument used for an instant model plot

if nargin < 2
    PLOT = 0;
end

[dy,dz,rho,type] = read_mackie2d_model(cfile);

if isempty(type)
    type = 'LINEAR';
end

ny = length(dy);
nzEarth = length(dz);
nzAir = 10;

dzAir(1:nzAir) = 0;
dzAir(1) = max(dz(1),10.0);
for k = 2:nzAir
   dzAir(k) = dzAir(k-1)*3;
end

grid.Dy = dy';
grid.Dz = [dzAir(end:-1:1) dz'];
grid.Ny = ny;
grid.Nz = nzEarth + nzAir;
grid.Nza = nzAir;

if findstr(type,'LOGE')
    cond.v = - rho;
else
    cond.v = 1./rho;
end
cond.AirCond = log(1e-10);
cond.paramType = type;
cond.grid = grid;

if PLOT
    [modelpath,modelname] = fileparts(cfile);
    str = ['2D Model: ' modelname];
    OPTIONS = struct('nYskip',0,'nZplot',nzEarth,'cax',[-4 0],'title',str);
    h = plotCond(cond,grid,OPTIONS);
    saveas(h,fullfile(modelpath,modelname),'png');
end