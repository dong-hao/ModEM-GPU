function [index,nx,ny,nz] = read_weerachai_control(fname)
% reads a 3D control file in Weerachai Siripunvaraporn's format
fid = fopen(fname);
line = fgetl(fid);
[n] = sscanf(line,'%d',[3 1]);
nx =   n(1);   ny  =   n(2);   nz  =   n(3);
k = 1;
index(1:nx,1:ny,1:nz) = -1;
while max(k) < nz
    tmp =   fgetl(fid);
    k   =   sscanf(tmp,'%d');
    if length(k)<2; k = [k k]; end
    tmp =   fscanf(fid,'%f',[ny nx]);
    for i = k(1):k(2)
       index(:,:,i) = tmp';
    end
    tmp =   fscanf(fid,'\n');
end
fclose(fid);