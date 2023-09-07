function status = write_mackie3d_model(fname,x,y,z,rho,nzAir,type,origin,rotation)
% writes a 3D resistivity model in Randie Mackie's format;
% allows for natural log resistivity by setting type = 'LOGE'
%  (c) Gary Egbert, 2007
%  open file for output
fid = fopen(fname,'w');
[nx, ny, nz] = size(rho);
if nargin <= 6
    type = ' ';
end
% output file
fprintf(fid, '%d %d %d %d %s %s\n', nx, ny, nz, nzAir, 'VALUES', type);
for j = 1:nx
    status = fprintf(fid,'%G ',x(j));
end
fprintf(fid, '\n');
for j = 1:ny
    status = fprintf(fid,'%G ',y(j));
end
fprintf(fid, '\n');
for j = 1:nz
    status = fprintf(fid,'%G ',z(j));
end
fprintf(fid, '\n');
for k = 1:nz
   fprintf(fid,'%d\n',k);
   for j = 1:ny
       % x index incremented fastest
       for i = 1: nx
        fprintf(fid,'%13.5E',rho(i,j,k));
       end    
       fprintf(fid, '\n');
    end
end
%  some WINGLINK garbage that we write out, then ignore ...
fprintf(fid, '%s \n', 'WINGLINK');
fprintf(fid, '%s \n', 'site');
fprintf(fid, '%d %d \n', 1,1);
%  add origin, rotation angle to end
if nargin <= 7
    origin = [0 0 0];
    rotation = 0;
end
fprintf(fid, '%d %d %d\n', origin/1000.0);
fprintf(fid, '%d\n', rotation);
status = fclose(fid);