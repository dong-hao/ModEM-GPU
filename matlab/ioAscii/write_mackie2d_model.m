function status = write_mackie2d_model(fname,y,z,rho,type)
% writes a 2D resistivity model in Randie Mackie's format
%  (c) Anna Kelbert, 2008
if nargin < 5
    type = ' ';
end
fid = fopen(fname,'w');
[ny, nz] = size(rho);
status = fprintf(fid,'%3d %3d %s\n',ny,nz,type);
for j = 1:ny
    status = fprintf(fid,'%11.4G',y(j));
    if (rem(j,10)==0) || (j==ny)
        fprintf(fid,'\n');
    end
end
for k = 1:nz
    status = fprintf(fid,'%11.4G',z(k));
    if (rem(k,10)==0) || (k==nz)
        fprintf(fid,'\n');
    end
end    
status = fprintf(fid,'1\n');
for k = 1:nz
    fprintf(fid,'  ');
    for j = 1:ny
        status = fprintf(fid,'%13.5E',rho(j,k));
        if (rem(j,19)==0) || (j==ny)
            fprintf(fid,'\n');
        end
    end
end
status = fclose(fid);