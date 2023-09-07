function status = write_meshtools3d_model(fname,Cond)
% writes a 3D conductivity model and mesh in the format
% suitable for MeshTools3D input

if strcmp(Cond.paramType,'LOGE')
    Cond.paramType = 'LINEAR';
    Cond.v = exp(Cond.v);
    Cond.AirCond = exp(Cond.AirCond);
end

nx = Cond.grid.Nx;
ny = Cond.grid.Ny;
nz = Cond.grid.NzEarth;
origin = 1000*Cond.grid.origin;
dx = Cond.grid.dx;
dy = Cond.grid.dy;
dz = Cond.grid.dz;
%origin(3) = sum(dz);

% mesh
fid = fopen([fname '.msh'],'w');
fprintf(fid, '%d %d %d\n', nx, ny, nz);
fprintf(fid, '%G %G %G\n\n', origin);
fprintf(fid,'%G %G %G %G %G\n',dx);
fprintf(fid,'\n\n');
fprintf(fid,'%G %G %G %G %G\n',dy);
fprintf(fid,'\n\n');
fprintf(fid,'%G %G %G %G %G\n',dz);
status = fclose(fid);

sigma = Cond.v;

%model
fid = fopen([fname '.con'],'w');
% for i = 1:nx 
%    for j = 1:ny
%        % z index incremented fastest
%        for k = 1:nz
%            fprintf(fid,'%13.5E\n',sigma(i,j,k));
%        end    
%        %fprintf(fid, '\n');
%     end
% end
% most likely, the above is the correct order...
% need to sort this out!
for j = 1:ny
   for i = 1:nx
       % z index incremented fastest
       for k = 1:nz
           fprintf(fid,'%13.5E\n',sigma(i,j,k));
       end    
       %fprintf(fid, '\n');
    end
end
status = fclose(fid);