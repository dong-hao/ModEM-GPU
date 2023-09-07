function newCond = interpCond_3D(newgrid,oldCond,bg)

% newCond = interpCond_3D(newgrid,oldCond,bg)
%
% Interpolates the conductivity parameter oldCond to the grid specified by
% newgrid. Argument bg gives the background absolute conductivity values to
% be used to replace NaN's (if any) - could be a scalar or a vector
% changing with depth. Use bg=0 if interpolating conductivity variations.

oldgrid = oldCond.grid;
oldv = oldCond.v;

oldx = oldgrid.origin(1)+cumsum([0; oldgrid.dx]);
oldy = oldgrid.origin(2)+cumsum([0; oldgrid.dy]);
oldz = oldgrid.origin(3)+cumsum([0; oldgrid.dz]);

newx = newgrid.origin(1)+cumsum([0; newgrid.dx]);
newy = newgrid.origin(2)+cumsum([0; newgrid.dy]);
newz = newgrid.origin(3)+cumsum([0; newgrid.dz]);

[oldY,oldX,oldZ] = meshgrid(oldy(1:end-1)+oldgrid.dy/2,...
    oldx(1:end-1)+oldgrid.dx/2,...
    oldz(1:end-1)+oldgrid.dz/2);

[newY,newX,newZ] = meshgrid(newy(1:end-1)+newgrid.dy/2,...
    newx(1:end-1)+newgrid.dx/2,...
    newz(1:end-1)+newgrid.dz/2);

newv = interp3(oldY,oldX,oldZ,oldv,newY,newX,newZ);

% replace NaN's, if any, with these values... use these for conductivity;
% use zero for conductivity perturbations
nz = length(newgrid.dz);
if isscalar(bg)
    bg(1:nz) = bg;
end
if strcmp(oldCond.paramType,'LOGE') && (min(bg) > 0)
    bg = log(bg);
end

for k = 1:nz
    temp = squeeze(newv(:,:,k));
    temp(isnan(temp)) = bg(k);
    newv(:,:,k) = temp;
end

newCond = oldCond;
newCond.grid = newgrid;
newCond.v = newv;