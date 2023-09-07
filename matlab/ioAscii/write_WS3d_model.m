function status = write_WS3d_model(fname,x,y,z,rho,nzAir,type,origin,rotation)
% writes a 3D resistivity model in Weerachai Siripunvaraporn's format;
% allows for natural log resistivity by setting type = 'LOGE'
%  (c) Anna Kelbert, 2009
%  open file for output
fid = fopen(fname,'w');
[nx, ny, nz] = size(rho);
if nargin <= 6
    type = ' ';
end
% output file
comment = 'Written by Matlab write_WS3d_model script'; 
fprintf(fid, '# %s\n', comment);
fprintf(fid, '%d %d %d %d %s\n', nx, ny, nz, 0, type);
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
   fprintf(fid,'\n');
   for j = 1:ny
       % x index incremented fastest
       for i = nx:-1:1
        fprintf(fid,'%15.5E',rho(i,j,k));
       end    
       fprintf(fid, '\n');
    end
end
%  add origin, rotation angle to end
if nargin <= 7
    origin = [-sum(x)/2 -sum(y)/2 0];
    rotation = 0;
end
fprintf(fid, '%d %d %d\n', origin);
fprintf(fid, '%d\n', rotation);
status = fclose(fid);