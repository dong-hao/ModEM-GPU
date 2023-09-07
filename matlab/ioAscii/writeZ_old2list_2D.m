fdir = '.';
filt = [fdir '/*.imp'];
[filename, pathname] = uigetfile(filt, 'Impedance File: Old Format');
fdata = [pathname filename];
[data,header,isign,units] = readZ_2D_old(fdata);
listname = [filename(1:end-4) '.dat'];
writeZ_2D(listname,data,header,isign,units);