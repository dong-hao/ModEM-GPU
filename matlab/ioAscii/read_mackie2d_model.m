function [y,z,rho,type] = read_mackie2d_model(fname)
% reads a 2D resistivity model in Randie Mackie's format;
% allows for natural log resistivity format
% assumes linear resistivity by default
fid = fopen(fname);
line = fgetl(fid);
[n] = sscanf(line,'%d %d',[2 1]);
if findstr(line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end
ny  =   n(1);   nz  =   n(2);
y   =   fscanf(fid,'%f',ny);
z   =   fscanf(fid,'%f',nz);
tmp =   fscanf(fid,'%d',1);
rho =   fscanf(fid,'%f',[ny nz]);
fclose(fid);