function status = write_weerachai_control(fname,index,format)
% writes a 3D control file in Weerachai Siripunvaraporn's format
% if format = 'x', x index is incremented fastest (like in our covariance);
% if format = 'y', y index is incremented fastest (like in WS code)
%  (c) Anna Kelbert, 2008
if nargin < 3
    format = 'y';
end
%  open file for output
fid = fopen(fname,'w');
[nx, ny, nz] = size(index);
% output file
fprintf(fid, ' %d %d %d\n', nx, ny, nz);
for k = 1:nz
   fprintf(fid,' %d %d\n',k,k);
   if strcmp(format,'y')
   for i = 1:nx
       % y index incremented fastest
       for j = 1:ny
        fprintf(fid,' %d',index(i,j,k));
       end    
       fprintf(fid, '\n');
   end
   elseif strcmp(format,'x')
   for j = 1:ny
       % x index incremented fastest
       for i = 1:nx
        fprintf(fid,' %d',index(i,j,k));
       end    
       fprintf(fid, '\n');
   end       
   end
end
status = fclose(fid);