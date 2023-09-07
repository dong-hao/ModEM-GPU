function cov = readCov_3D(cfile)

%  Usage:  cov = readCov_3D(cfile)
%  
%  Reads a 3-D covariance matrix from the Mod3DMT format.
%  In the covariance matrix, air = 0 and ocean = 9.
%
%  (c) Anna Kelbert, 2009

fid = fopen(cfile,'r');

% skip the header
for i = 1:16
    fgetl(fid);
end

% read the model size
nx = fscanf(fid,'%d',1);
ny = fscanf(fid,'%d',1);
nz = fscanf(fid,'%d',1);
disp(['Model size: ' num2str(nx) 'x' num2str(ny) 'x' num2str(nz)]);

% initialize the covariance
cov = zeros(nx,ny,nz);

% skip the smoothing parameters
fscanf(fid,'%G',nz);
fscanf(fid,'%G',nz);
fscanf(fid,'%G',1);
fscanf(fid,'%d',1);

% skip the rules
nrules = fscanf(fid,'%d',1);
for i = 1:nrules
    fgetl(fid);
end
fscanf(fid, '\n');

% read the covariance
k = [0 0];
while max(k) < nz
    tmp = fgetl(fid);
    k   = sscanf(tmp,'%d');
    tmp = fscanf(fid,'%d',[ny nx]);
    for K = k(1):k(2)
       cov(:,:,K) = flipud(tmp');
    end
    fscanf(fid,'\n');
end

% close file
fclose(fid);
