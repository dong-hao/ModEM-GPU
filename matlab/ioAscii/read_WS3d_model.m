function [x,y,z,rho,nzAir,type,origin,rotation] = read_WS3d_model(fname)
% reads a 3D resistivity model in Weerachai Siripunvaraporn's "0" format
fid = fopen(fname);
line = fgetl(fid); % comment
line = fgetl(fid);
[n] = sscanf(line,'%d',[4 1]);
if findstr(line,'LOGE')
    type = 'LOGE';
else
    type = 'LINEAR';
end
nx =   n(1);   ny  =   n(2);   nz  =   n(3);
nzAir = 6;
x   =   fscanf(fid,'%f',nx);
y   =   fscanf(fid,'%f',ny);
z   =   fscanf(fid,'%f',nz);
k = 1;
rho(1:nx,1:ny,1:nz) = 0;
while k <= nz
    for j = 1:ny
        for i = nx:-1:1
            rho(i,j,k) = fscanf(fid,'%f',1);
        end
    end
    k = k+1;
end
line = fgetl(fid); % read to the end of line
origin = [-sum(x)/2 -sum(y)/2 0];
rotation = 0;
while 1
    line = fgetl(fid);
    if ~ischar(line); 
        fclose(fid);
        return; 
    end
    [n] = sscanf(line,'%f',[3 1]);
    if length(n)==3
        origin = n;
    else
        [n] = sscanf(line,'%f',1);
        if length(n)==1
            rotation = n; 
        end
    end
end