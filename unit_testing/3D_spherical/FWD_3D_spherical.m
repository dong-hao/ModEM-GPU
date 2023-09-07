%%   load in model parameter, and grid
%   using mTrue.cpr as input, results can be compared to e0.soln
%   This can be obtained by running 
%     Mod3DMT -F -F mTrue.cpr Template_mTrue.dat out.dat e0.soln
%   For testing on other models similar procedures can be used, though
%   conceivably the data template file would have to be modified.
modelFile = 'CylinderModel+QuarterSpace+Spherical.rho';
%modelFile = 'cascadia.ws';
%modelFile = 'commemi3d2.ws';
%MP = 'SG';
MP = 'GDE Spherical';
%modelFile = 'mTrueS.cpr';
%   spherical test case
%modelFile = 'model_true2.mod';
%modelFile = 'model_true_large_4_fineCenter.mod';
%modelFile = 'mNew_4_48.cpr';%
%modelFile = 'Unif6x6x7.cpr';
%nAirLayers = 5;
%airTop = 1e+5;
%modelFile = 'mNew_4.cpr';
nAirLayers = 12;
airTop = 1e+6;
%modelFile = 'mUnif_4.cpr';
%modelFile = 'm_100_24.cpr';
%modelFile = 'm_100_16unif.cpr';
%MP = 'GDE Spherical';
%Solver = 'Direct';
%Solver = 'Iter';
Solver = 'Iter2';
%Solver = 'IterDiv';
%Solver = 'Will';
uniform=false;

%   this is all of the periods for MT3D workshop
%T = [.01,.018,0.032,0.056,.1,.18,0.32,.56,1,1.8,3.2,5.6,10,18,32,56,100,180,320,560,...
%    1000,1800,3200,5600,10000];
%   here is a subset
%T = [.01,.03,.1,.3,1,3.2,10,32,100];%,320,1000,3200,10000];
dat = mtdata.read('CylinderModel_6per_Spherical_Template.dat','list','[mV/km]/[nT]','Full_Impedance');
T = dat.periods;
nPer = length(T);
location = dat.d{1}.siteLoc; % site position
nMode = 2;
XY = {'X','Y'};
SOLNS = struct('grid',[],'Periods',T,'Modes',{XY},'E',[],'S',[],'Diagnostics',[]);
SOLNS.E = cell(nMode,nPer);
%   this block is specific to a particular model parameter implementation
switch MP  
    case 'SG'
        ModPar = readVec(TModelParameterCell3D_EC_SG(),modelFile);
        %  to make this clear, I am creating a separate grid object using the
        %   grid from the model parameterization, m.ParamGrid, which now
        %   does not include air layers
        ModPar.grid = setAirLayers(ModPar.ParamGrid,'maxheight',airTop,...
            'nlayers',nAirLayers);
           % 'method','mirror','nlayers',10);
        ModOp = TModelOperator3DC_SG(ModPar.grid);
    case 'GDE Spherical'
        ModPar = readVec(TModelParameterCell3D_ES_SG(),modelFile);
        %  to make this clear, I am creating a separate grid object using the
        %   grid from the model parameterization, m.ParamGrid
        ModPar.grid = setAirLayers(ModPar.ParamGrid,'maxheight',1e+6,...
            'nlayers',12);%,'method','mirror','nlayers',10);
        ModOp = TModelOperator3DS_SG(ModPar.grid);
end
%ModPar.AirCond = log(1e-10);
SOLNS.grid = ModPar.grid;
SOLNS.grid.rotation = 0;
SOLNS.S = ModPar;
if uniform
    ModPar.v = ModPar.v(1,1,1)*ones(size(ModPar.v)); %#ok<UNRCH>
end
%%
%    create forward modeling operator: probably in normal usage just do grid
%    dependent setup; setup for model parameter and frequency separately
switch Solver
    case 'Direct'
        fwd = TForwardModel(ModOp,ModPar);            % direct solver
    case 'Iter'
        fwd = TForwardModel_IT(ModOp,ModPar);       % iterative solver without div correction
    case 'Iter2'
        fwd = TForwardModel_IT2(ModOp,ModPar);  % iterative solver of modified
        %   system: grad div added in air,
        %   sig grad div sig in  Earth
    case 'IterDiv'
        fwd = TForwardModel_ITdiv(ModOp,ModPar);  % iterative solver with div correction
    case 'Will'
        fwd = TForwardModel_WillGDE(ModOp,ModPar);  % iterative solver with 
end

fwd.setEquations;

if isa(fwd,'TForwardModel_IT')
   SOLNS.Diagnostics = cell(nMode,nPer);
   fwd.setSolverOptions('EMSolverTolerance',1e-10);
   fwd.setSolverOptions('EMSolverMethod','BCG') %   opitons: BCG, QMR, GMRES
   fwd.setSolverOptions('EMSolverPrecond','ILU0') %   opitons: ILU0,  GS, SSOR
   fwd.setSolverOptions('Initialization','none')%   options: none, 1D, 2DTE, E0
     %   somewhat astoundingnly "none" seems slightly faster than the other
     %   two options!
     %   (probably does not make any difference to convergence, just saves
     %   the   calculations of starting model)

   fwd.setSolverOptions('IterPerDivCor',80)
   fwd.setSolverOptions('MaxDivCor',30)
   fwd.setSolverOptions('DivCorTolerance',1e-7)
   fwd.setSolverOptions('DivCorPreCond','ICD')
end

if isa(fwd,'TForwardModel_IT2')
    fwd.nodeORedge = 'node';
    fwd.variable = true;
end
%%  ultimately make a loop over periods and modes
E = TVector3D_SG(ModPar.grid);
tic
 for iperiod = 1:nPer
     fprintf('%s %f\n','Period',T(iperiod))
     omega = 2*pi/T(iperiod);   %  angular frequency
     if iperiod== 1
         fwd.setFreqCond(omega,ModPar);
     else
         %    only need to set frequency, not model parameter
         fwd.setFreqCond(omega)
     end
     if isa(fwd,'TForwardModel_IT')
         fwd.setPreConditioner;
     end
     for imode = 1:nMode
         polarization = XY{imode};
         fprintf('%s %s\n','Polarization',polarization);
         %   set up source term
         src = TSourceMT_1DBC(ModPar);
         src.SetSourceParams(omega,polarization);
         %   solve
         [e,diagnostics] = fwd.Solve(src);
         SOLNS.E{imode,iperiod} =  E.SetVector(e);
         SOLNS.Diagnostics{imode,iperiod} = diagnostics;
     end
   %  data = TMTsite(fwd,location,omega);
   %  Z = data.Impedance(SOLNS.E{1,1},SOLNS.E{2,1});
 end
 toc
 eval(['save FWD_SOLNS_6_' Solver '.mat SOLNS fwd location'])
 plotEMsoln(SOLNS)

 %% read ModEM outputs
dir = 'OUTPUT/2020-09-03-stable-SP2/';
%dir = 'OUTPUT/2020-09-03-gridelements/';
dat = mtdata.read([dir 'CylinderModel_6per.dat'],'list','[mV/km]/[nT]','Full_Impedance');
T = dat.periods;
nPer = length(T);

% quick compare
isite = 2000; iperiod = 6; 
disp(['Site: ' char(dat.d{1}.siteChar(isite)) ' ' ...
    num2str(dat.d{1}.lat(isite)) ' ' num2str(dat.d{1}.lon(isite))]);
pred = dat.zero;
site = TMTsite(fwd,location(isite,:));
func = TDataFunc(site);
[Z,H] = func.ComputeData(SOLNS.E{1,iperiod},SOLNS.E{2,iperiod});
factor = ImpUnits('[Ohm]','[mV/km]/[nT]');
pred.d{iperiod} = pred.d{iperiod}.set(isite,Z);
disp('FORTRAN: '); dat.d{iperiod}.TF(isite,:)
disp(' MATLAB: '); conj(pred.d{iperiod}.TF(isite,:))*factor*(1000/(4^(iperiod-1)))

datplot = mtdataplot(dat);
datplot.plot('ZYX',1);

% load the E-fields
SOLNSFOR = readE([dir 'CylinderModel_6per.soln'],6);
plotEMsoln(SOLNSFOR);

%ExyzPlot;
SOLNSgridelements = readE('../3D/CylinderModel_6per_gridelements.soln',6);
plotEMsoln(SOLNSgridelements)

dir = 'OUTPUT/2020-09-03-gridelements/';
SOLNSgridelementsS = readE([dir 'CylinderModel_6per.soln'],6);
plotEMsoln(SOLNSgridelementsS)

diff = SOLNSstable.E{1,6}-SOLNSgridelements.E{1,6};


 %% compute and save the impedances in a list
pred = dat;
for isite = 1:length(location)
    for iperiod = 1:nPer
        site = TMTsite(fwd,location(isite,:));
        func = TDataFunc(site);
        [Z,H] = func.ComputeData(SOLNS.E{1,iperiod},SOLNS.E{2,iperiod});
        pred.d{iperiod} = pred.d{iperiod}.set(isite,Z);
    end
end
pred.write('CylinderModel_6per_Matlab.dat','list');

%% add noise and save again
pred = mtdata.read('CylinderModel_1per.dat','list','[mV/km]/[nT]','Full_Impedance');
pred_noisy = pred.addNoise(0.05);
pred_noisy.write('CylinderModel_1per_noisy.dat','list');

pred = mtdata.read('CylinderModel_6per.dat','list','[mV/km]/[nT]','Full_Impedance');
pred_noisy = pred.addNoise(0.05);
pred_noisy.write('CylinderModel_6per_noisy.dat','list');

%% add noise and save again
pred_noisy = pred.addNoise(0.05);
pred_noisy.write('CylinderModel_6per_Matlab_noisy.dat','list');

%% compute and save the impedances - use individual TFs - takes time
% tflist = dat.mtdata2tflist;
% for isite = 1:length(location)
%     Z = zeros(2,2,nPer);
%     for iperiod = 1:nPer
%         site = TMTsite(fwd,location(isite,:));
%         func = TDataFunc(site);
%         [Z(:,:,iperiod),H] = func.ComputeData(SOLNS.E{1,iperiod},SOLNS.E{2,iperiod});
%     end
%     tflist{isite} = tflist{isite}.set(T,Z);
% end
% pred = mtdata.tflist2mtdata(tflist,'Full_Impedance',T);
% pred.write('FWD_3D_spherical.dat','list');