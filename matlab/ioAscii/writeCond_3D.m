function [status] = writeCond_3D(cfile,cond,format)
%
% Usage:  [status] = writeCond_3D(cfile,cond,format)
%
%  writes ModelParam object, provided as structure
%  status is total number of bytes written
%
%  format = 1 means Mackie's format (default)
%  format = 2 means WS's format

if nargin < 3
    format = 1;
end

nzAir = cond.grid.NzAir;
dx = cond.grid.dx;
dy = cond.grid.dy;
dz = cond.grid.dz;

if isfield(cond.grid,'origin')
    origin = cond.grid.origin;
else
    origin = [0 0 0];
end

if isfield(cond.grid,'units')
    if strcmp(cond.grid.units,'km')
        % convert everything to meters!
        dx = 1000*dx;
        dy = 1000*dy;
        dz = 1000*dz;
        origin = 1000*origin;
    end
end        

if isfield(cond.grid,'rotation')
    rotation = cond.grid.rotation;
else
    rotation = 0;
end

if findstr(cond.paramType,'LOGE')
    type = 'LOGE';
    rho = - cond.v;
elseif findstr(cond.paramType,'LOG10')
    type = 'LOG10';
    rho = - cond.v;
else
    type = 'LINEAR';
    rho = 1./(cond.v);
end

if format == 1
    status = write_mackie3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
elseif format == 2
    status = write_WS3d_model(cfile,dx,dy,dz,rho,nzAir,type,origin,rotation);
else
    error('Unknown format');
end
