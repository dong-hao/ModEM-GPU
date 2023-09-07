prior = 1000;
Cond = readCond_3D('Yellowstone_10km_pr.ws',2);
Cond.grid.Ny = 66;
Cond.grid.xy = 2*Cond.grid.dx;
Cond.grid.dy = 2*Cond.grid.dy;
Cond.grid.dy = [Cond.grid.dy(1:5); Cond.grid.dy(20:end-19); Cond.grid.dy(end-4:end)];
Cond.v = (1/prior)*ones(Cond.grid.Nx,Cond.grid.Ny,Cond.grid.NzEarth);
Cond.grid.origin(1) = - mean(cumsum(Cond.grid.dx)); 
Cond.grid.origin(2) = - mean(cumsum(Cond.grid.dy)); 
[status] = writeCond_3D(['rubiks_20km_' num2str(prior) 'ohmm_pr.ws'],Cond,2);

cube_bg = 100;
a1 = 0.1;
a3 = 0.001;
cube1 = 7:22;
cube2 = 26:41;
cube3 = 45:60;
Cond0 = Cond;
Cond0.v = (1/cube_bg)*ones(Cond.grid.Nx,Cond.grid.Ny,Cond.grid.NzEarth);
% layer 1
Cond0.v(cube1,cube1,1:16) = a3;
Cond0.v(cube2,cube1,1:16) = a1;
Cond0.v(cube3,cube1,1:16) = a3;
Cond0.v(cube1,cube2,1:16) = a1;
Cond0.v(cube2,cube2,1:16) = a3;
Cond0.v(cube3,cube2,1:16) = a1;
Cond0.v(cube1,cube3,1:16) = a3;
Cond0.v(cube2,cube3,1:16) = a1;
Cond0.v(cube3,cube3,1:16) = a3;
% layer 2
Cond0.v(cube1,cube1,21:28) = a1;
Cond0.v(cube2,cube1,21:28) = a3;
Cond0.v(cube3,cube1,21:28) = a1;
Cond0.v(cube1,cube2,21:28) = a3;
Cond0.v(cube2,cube2,21:28) = a1;
Cond0.v(cube3,cube2,21:28) = a3;
Cond0.v(cube1,cube3,21:28) = a1;
Cond0.v(cube2,cube3,21:28) = a3;
Cond0.v(cube3,cube3,21:28) = a1;
% layer 3
Cond0.v(cube1,cube1,30:35) = a3;
Cond0.v(cube2,cube1,30:35) = a1;
Cond0.v(cube3,cube1,30:35) = a3;
Cond0.v(cube1,cube2,30:35) = a1;
Cond0.v(cube2,cube2,30:35) = a3;
Cond0.v(cube3,cube2,30:35) = a1;
Cond0.v(cube1,cube3,30:35) = a3;
Cond0.v(cube2,cube3,30:35) = a1;
Cond0.v(cube3,cube3,30:35) = a3;
[status] = writeCond_3D('rubiks_20km.ws',Cond0,2);

x = cumsum(Cond.grid.dx');
y = cumsum(Cond.grid.dy');
z = cumsum(Cond.grid.dz');
x = [0 x];
y = [0 y];
z = [0 z];
xctr = x(1:end-1) + diff(x)/2 + Cond.grid.origin(1);
yctr = y(1:end-1) + diff(y)/2 + Cond.grid.origin(2);
xcnr = x(1:end-1) + Cond.grid.origin(1);
ycnr = y(1:end-1) + Cond.grid.origin(2);

% create a data point at all cell centers except padding
xpad = 5;
ypad = 5;
xdiff = 4; % create sites every xdiff cells (latitude)
ydiff = 4; % create sites every ydiff cells (longitude)
siteLoc = [];
siteChar = '';
for i = xpad+1:xdiff:Cond.grid.Nx-xpad+1
    for j = ypad+1:ydiff:Cond.grid.Ny-ypad+1
        siteLoc = [siteLoc; xcnr(i) ycnr(j) 0.0];
        siteChar = [siteChar; sprintf('%03d',i) '-' sprintf('%03d',j)];
    end
end

nSites = size(siteLoc,1);

%   PERIODS:
T1 = 10;
T2 = 10000;
nPer = 12;
if nPer > 1
    dt = (log10(T2/T1))/(nPer-1);
    periods = 10.^[log10(T1):dt:log10(T2)];
else
    periods = T1;
end

nImp = 6;
compChar = ['Re(Zxx)';'Im(Zxx)';'Re(Zxy)';'Im(Zxy)'; ...
        'Re(Zyx)';'Im(Zyx)';'Re(Zyy)';'Im(Zyy)'; ...
        'Re(Tx) ';'Im(Tx) ';'Re(Ty) ';'Im(Ty) '];

Z = zeros(nSites,nImp)+1i*zeros(nSites,nImp);
Zerr = ones(nSites,nImp);

allData = cell(nPer,1);
for j=1:nPer
    allData{j} = struct('T',periods(j),...
        'nComp',nImp*2,...
        'compChar',compChar,...
        'siteLoc',siteLoc*1000,...
        'siteChar',siteChar,...
        'Z',Z,'Zerr',Zerr,...
        'Cmplx',1);
end
%   end data vector template setup
[status] = writeZ_3D('rubiks_Template.dat',allData);

% NOW COMPUTE RESPONSES EXTERNALLY; 
% mpirun -n 13 ../Mod3DMT_MPI -F rubiks_20km.ws rubiks_Template.dat rubiks_20km_fwd.dat
% then load and add noise
FracError = 0.03;
[data] = readZ_3D('rubiks_20km_fwd.dat');
display(['adding noise ... FracError = ' num2str(FracError)]) 
data = addNoise_3D(data,FracError);
[status] = writeZ_3D('rubiks_20km_3%.dat',data);
