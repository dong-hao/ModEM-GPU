fdir = '.';
filt = [fdir '/*.imp'];
[filename, pathname] = uigetfile(filt, 'Impedance File: Old Format');
fdata = [pathname filename];
[data,header,isign,units] = readZ_3D_old(fdata);
SYNTHETIC = 0;
if ~ SYNTHETIC %  if GG locations are not supplied, create them
    for i=1:length(data)
        [data{i}.lat,data{i}.lon] = ...
            xy2latlon(data{i}.siteLoc(:,1),data{i}.siteLoc(:,2), ...
                      data{i}.origin(1),data{i}.origin(2),'m');
    end
end        
listname = [filename(1:end-4) '.dat'];
writeZ_3D(listname,data,header,isign,units);