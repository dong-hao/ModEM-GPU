function [x,y,z,rho,nzAir,type,origin,rotation] = read_mackie3d_model(fname,block)
% reads a 3D resistivity model in Randie Mackie's format
% if 'block' is specified and true, no integer indices are read
% and the model is transposed (provides support for WS's models)
if nargin < 2
    block = 0;
end
fid = fopen(fname);
line = fgetl(fid);
[n] = sscanf(line,'%d',[4 1]);
if findstr(line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end
nx =   n(1);   ny  =   n(2);   nz  =   n(3);
nzAir = n(4);
x   =   fscanf(fid,'%f',nx);
y   =   fscanf(fid,'%f',ny);
z   =   fscanf(fid,'%f',nz);
k = 1;
rho(1:nx,1:ny,1:nz) = 0;
if block
    for i = 1:nz
        tmp = fscanf(fid,'%f',[nx ny]);
        rho(:,:,i) = flipud(tmp);
    end
else
while max(k) < nz
    tmp =   fscanf(fid,'%s\n',1);
    k   =   sscanf(tmp,'%d');
    if length(k)<2; k = [k k]; end
    tmp =   fscanf(fid,'%f',[nx ny]);
    for i = k(1):k(2)
       rho(:,:,i) = tmp;
    end
end
end
origin = [0 0 0];
rotation = 0;
line = fgetl(fid); if ~ischar(line); return; end
line = fgetl(fid); if ~ischar(line); return; end
line = fgetl(fid); if ~ischar(line); return; end
line = fgetl(fid); if ~ischar(line); return; end
line = fgetl(fid); if ~ischar(line); return; end
[n] = sscanf(line,'%f',[3 1]);
if length(n)==3; origin = 1000.0*n; end
line = fgetl(fid); if ~ischar(line); return; end
[n] = sscanf(line,'%f',1);
if length(n)==1; rotation = n; end
fclose(fid);