filt = './*.ws; *.cpr; *.rho; *.prm';
[filename, pathname] = uigetfile(filt, 'Resistivity Model File');

%  read in chosen file
cfile = [pathname filename];

Cond = readCond_3D(cfile,2);

% convert to log10
if strcmp(Cond.paramType,'LOGE')
    Cond.paramType = 'LOG10';
    Cond.v = Cond.v / log(10);
    Cond.AirCond = Cond.AirCond / log(10);
elseif strcmp(Cond.paramType,'LINEAR')
    Cond.paramType = 'LOG10';
    Cond.v = log10(Cond.v);
    Cond.AirCond = log10(Cond.AirCond);
end

filt = './*.dat; *.res; *.imp';
[filename, pathname] = uigetfile(filt, 'Load Data File');
Zfile = [pathname filename];

[AllData] = readZ_3D(Zfile,'[mV/km]/[nT]');

plot3Dview
caxis([-3.5 -0.5]);

print('-dpng','-r300', cfile(1:end-3));
