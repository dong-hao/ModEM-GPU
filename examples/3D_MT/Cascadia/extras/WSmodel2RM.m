% convert from Weerachai control file format to ModEM covariance

% read Weerachai control file
cfile = 'original/control_new';
ctl = read_weerachai_control(cfile); % 80 x 78 x 34

% find ocean and update index
ctl(ctl==1) = 9;

% find the rest of the model and update index
ctl(ctl==0) = 1;

% for i=1:size(ctl,3)
%     ctl(:,:,i) = flipud(ctl(:,:,i));
% end

% rearrange array to fit nx values in a row format
% ctl = permute(ctl,[2 1 3]);

% write file
cfile = 'temp.cov';
write_weerachai_control(cfile,ctl);

ctlhalf = ctl([1:8 9:2:72 73:80],[1:7 8:2:71 72:78],:); % 48 x 46 x 34

cfile = 'temp_half.cov';
write_weerachai_control(cfile,ctlhalf);


% read Weerachai prior model file
cfile = 'original/cascad_half_model.00';
Cond = readCond_3D(cfile,2);

% update the model for ModEM
Cond.paramType = 'LOGE';
Cond.v = log(Cond.v);
Cond.AirCond = log(Cond.AirCond);
Cond.grid.NzAir = 6;
Cond.grid.origin = [-728840 -795244 0];
Cond.grid.rotation = 0;

% write model in Mackie's format
cfile = 'cascad_prior.rm';
writeCond_3D(cfile,Cond,1);

% make a half size model
CondHalf = Cond;
dx = Cond.grid.dx;
dy = Cond.grid.dy;
CondHalf.grid.Nx = 48;
CondHalf.grid.dx = [dx(1:7); 18000; 18000; 24000*ones(1,30)'; 18000; 18000; dx(74:80)];
CondHalf.grid.Ny = 46;
CondHalf.grid.dy = [dy(1:7); 24000*ones(1,32)'; dy(72:78)];
CondHalf.grid.origin = [-728840 -795244 0];
CondHalf.grid.rotation = 0;
CondHalf.v = Cond.v([1:8 9:2:72 73:80],[1:7 8:2:71 72:78],:); % 48 x 46 x 34

% write half size model in WS's format
cfile = 'cascad_half_prior.ws';
writeCond_3D(cfile,CondHalf,2);

% and in Mackie's format
cfile = 'cascad_half_prior.rm';
writeCond_3D(cfile,CondHalf,1);

% plot the linear conductivity at any time using the lines below
options.slice = 'Z';
options.Np = 1;
[Nx,Ny,Nz] = size(Cond.v);
options.iXlim(1) = 1;
options.iXlim(2) = Nx+1;
options.iYlim(1) = 1;
options.iYlim(2) = Ny+1;
options.iZlim(1) = 1;
options.iZlim(2) = Nz+1;
CondPlotSet(exp(Cond.v),Cond.grid,options);