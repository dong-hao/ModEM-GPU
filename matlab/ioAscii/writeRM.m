function [] = writeRM(cfile, dx, dy, dz, nzAir,rho, origin, rotation)
% usage: writeRM(cfile, dx, dy, dz, rho, origin, rotation);
%  Does not now support writing out in resistivity code format ...

%  open file for output
fid = fopen(cfile,'w');

[nx, ny, nz] = size(rho);

% output file
fprintf(fid, '%d %d %d %d %s\n', nx, ny, nz, nzAir, 'VALUES');
i = 1;
while i <= nx
    fprintf(fid, '%d ', dx(i));
    i = i + 1;
end
fprintf(fid, '\n');
i = 1;
while i <= ny
    fprintf(fid, '%d ', dy(i));
    i = i + 1;
end
fprintf(fid, '\n');
i = 1;
while i <= nz
    fprintf(fid, '%d ', dz(i));
    i = i + 1;
end
fprintf(fid, '\n');
for k = 1:nz
   fprintf(fid,'%d\n',k);
   for j = 1:ny
       % x index incremented fastest
       for i = 1: nx
        fprintf(fid,'%d ',rho(i,j,k));
       end    
       fprintf(fid, '\n');
    end
end
%  some WINGLINK garbage that we write out, then ignore ...
fprintf(fid, '%s \n', 'WINGLINK');
fprintf(fid, '%s \n', 'site');
fprintf(fid, '%d %d \n', 1,1);
%  add origin, rotation angle to end
fprintf(fid, '%d %d %d\n', origin);
fprintf(fid, '%d\n', rotation);
status = fclose('all');
