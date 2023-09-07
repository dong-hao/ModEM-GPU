
%%
dx = 10*ones(40,1);
dy = 10*ones(60,1);
dz = xygrid.logz(0,100,36,1.2);
dz2 = xygrid.logz(100,1000,8,1.5);
mygrid = xygrid(dx,dy,[dz;dz2]);

npad = 4; 
mygrid = mygrid.pad('NSEW',npad);

ind = mygrid.mask('cylinder',[5 15],[0 100],100,5);

sigma = xymodel(mygrid);
sigma.v = -2*ones(mygrid.nx,mygrid.ny,mygrid.nzEarth);
sigma = log10(sigma);
sigma.v(:,1:20+npad,:) = -1;
sigma.v(ind==5) = 1;
sigma.uiplot
sigma = loge(sigma);
sigma.write('CylinderModel+QuarterSpace.rho','WS');

prior = sigma;
prior.v(:,:,:) = log(0.01);
prior.write('Halfspace.rho','WS');

% projection copied from Paul's MCR, all corrected
utmstruct = defaultm('tranmerc');
utmstruct.origin = [45 -90 0];
utmstruct.scalefactor = 0.9996;
utmstruct.geoid = almanac('earth','wgs84','meters');
utmstruct.falseeasting=500000;
mstruct = defaultm(utmstruct);

% create data locations at nodes
per = 10.^(0.6:0.6:4);
dat = mtdata('Full_Impedance',per,'xy');
dat = dat.gridNodes(mygrid,per,'all',mstruct);
dat = zero(dat);
dat.header = 'Cylinder Model - Data at Nodes';
dat.write('CylinderModel_6per_Template.dat','list');

% ok, make also a single period version for quicker testing...
per = 1000;
dat = mtdata('Full_Impedance',per,'xy');
dat = dat.gridNodes(mygrid,per,'all',mstruct);
dat = zero(dat);
dat.header = 'Cylinder Model - Data at Nodes';
dat.write('CylinderModel_1per_Template.dat','list');

% not using a correct projection - just transforming the grid!
latgrid = sigma.grid.llgrid(mstruct);

% create an identical spherical coords model, no conversion
llsigma = llmodel(latgrid);
llsigma = loge(llsigma);
llsigma.v = sigma.v;
llsigma.uiplot;
llsigma = llsigma.setLimits;
llsigma.write('CylinderModel+QuarterSpace+Spherical.rho','ModEM');

llprior = llsigma;
llprior.v(:,:,:) = log(0.01);
llprior.write('Halfspace+Spherical.rho','ModEM');

% and using that new grid to create a new data set
per = 10.^(0.6:0.6:4);
dat = mtdata('Full_Impedance',per,'ll');
dat = dat.gridNodes(latgrid,per,'all');
dat = zero(dat);
dat.header = 'Cylinder Model - Data at Nodes';
dat.write('CylinderModel_6per_Spherical_Template.dat','list');

% and make also a single period version for quicker testing...
per = 1000;
dat = mtdata('Full_Impedance',per,'ll');
dat = dat.gridNodes(latgrid,per,'all');
dat = zero(dat);
dat.header = 'Cylinder Model - Data at Nodes';
dat.write('CylinderModel_1per_Spherical_Template.dat','list');

% also create a "correctly converted" spherical coords model, BAD!
llsigma0 = sigma.llmodel(mstruct);
llsigma0 = llsigma0.setLimits;
llsigma0 = llsigma0.fillNaNs;
llsigma0.uiplot;

modelname = 'CylinderModel+QuarterSpace';
sigma.plot('depth',10); delta=1e-6; caxis([-4-delta 0+delta]);
print('-djpeg','-r300',[modelname '_depth_' num2str(10) 'km.jpg']);

modelname = 'CylinderModel+QuarterSpace+Spherical';
llsigma.plot('depth',10); delta=1e-6; caxis([-4-delta 0+delta]);
print('-djpeg','-r300',[modelname '_depth_' num2str(10) 'km.jpg']);

%% TEST #1
% ../stable/f90/Mod3DMT -F CylinderModel+QuarterSpace.rho CylinderModel_6per_Template.dat CylinderModel_6per.dat CylinderModel_6per.esoln