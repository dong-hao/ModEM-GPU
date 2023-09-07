function [allData,info,units,isign] = readZ_2D(cfile)
%  Usage: [allData,info,units,isign] = readZ_2D(cfile)
%   info, units and isign are optional output arguments.
%   This information could be part of allData structure.
%   reads data from ascii impedance file, returns 
%    as cell array allData
%   There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   NB  Assuming complex data
%  (c) Anna Kelbert, 2008

fid = fopen(cfile,'r');

%  first read the header
record = textscan(fid,'%*[^:]%*c%[^\n]',3);
info = record{1}{1};
units = record{1}{2};
isign = str2num(record{1}{3});

%  record # 1: number of transmitters
nTx = fscanf(fid,'%d\n',1);

for j = 1:nTx
%   for each period/mode ....

   %  record # 1 in transmitter block: period, mode, number of sites
   record = textscan(fid,'%f %2c %d\n',1);
   T = record{1};
   MODE = record{2};
   nSites = double(record{3});
   nComp = 2;

   %  record # 2 in transmitter block: site locations
   record = fscanf(fid,'%f',[nSites,2]);
   siteLoc = record;
   record = fscanf(fid,'\n');

   %  read and ignore the header of the data block
   header = fgetl(fid);

   %  records # 3 & #4 in transmitter block: data and error bars
   siteChar = '';
   for k = 1:nSites
       idx = 1:2:nComp-1;
       record = fscanf(fid,'%s',1);
       siteChar = strvcat(siteChar, record);
       record = fscanf(fid,'%f',nComp);
       Z(k,:) = record(idx)+i*record(idx+1);
       record = fscanf(fid,'%f',nComp);
       Zerr(k,:) = record(idx);
   end

   allData{j} = struct('T',T,'Cmplx',1,...
	'Mode',MODE,'siteLoc',siteLoc,'siteChar',siteChar,'Z',Z,'Zerr',Zerr);

   clear siteChar
end
fclose(fid);
