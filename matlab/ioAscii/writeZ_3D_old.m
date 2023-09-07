function [status] = writeZ_3D(cfile,allData,info,units,isign)
%  Usage:  [status] = writeZ_3D(cfile,allData,info,units,isign);
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%  (c) Anna Kelbert, 2008

nTx = length(allData);
fid = fopen(cfile,'w');

%  description <= 80 char in length
if nargin < 3
    info = 'Synthetic 3D MT data written in Matlab';
end

%  assume SI units by default; alternative is [mV/km]/[nT] 
if nargin < 4
    units = 'Ohm';
end

%  sign convention is -1 by default
if nargin < 5
    isign = -1;
end

%  header: three lines followed by an empty line
status = fprintf(fid,'Description: %s\n',info);
status = fprintf(fid,'Units: %s\n',units);
status = fprintf(fid,'Sign convention: %d\n\n',isign);

%  record # 1: number of transmitters
status = fprintf(fid,'%d\n',nTx);

for j = 1:nTx
%   for each period/mode ....

   T = allData{j}.T;
   nComp = allData{j}.nComp;
   nSites = size(allData{j}.siteLoc,1);
   
   %  remove NaNs
   badSites = [];
   for k = 1:nSites
       if sum(isnan(allData{j}.Z(k,:)))>0
           badSites = [badSites k];
       end
   end
   goodSites = setdiff(1:nSites,badSites);

   %  record # 1 in transmitter block: period, mode, number of sites
   status = fprintf(fid,'%12.6E %5d %5d\n',T,nComp,length(goodSites));

   %  record # 2 in transmitter block: site locations
   siteLoc = allData{j}.siteLoc;
   for i = 1:3
       for k = goodSites
           status = fprintf(fid,'%12.3f ',siteLoc(k,i));
       end
       status = fprintf(fid,'\n');
   end

   %  header of the data & error block
   for i = 1:2:nComp-1
       if isfield(allData{j},'compChar')
         status = fprintf(fid,'%10s %s',' ',allData{j}.compChar(i,:));
         status = fprintf(fid,'%10s %s',' ',allData{j}.compChar(i+1,:));
       else
         status = fprintf(fid,'%12s Re',' ');
         status = fprintf(fid,'%12s Im',' ');
       end
   end
   status = fprintf(fid,'\n');

   %  records # 3 & #4 in transmitter block: data and error bars
   for k = goodSites
       if isfield(allData{j},'siteChar')
           siteChar = allData{j}.siteChar(k,:);
       else
           siteChar = sprintf('%03d',k);
       end
       status = fprintf(fid,'%10s',siteChar);
       Zri = allData{j}.Z;
       if allData{j}.Cmplx
           for i = 1:nComp/2
               status = fprintf(fid,'%15.6E %15.6E',real(Zri(k,i)),imag(Zri(k,i)));
           end
       else
           for i = 1:nComp
               status = fprintf(fid,'%15.6E',Zri(k,i));
           end
       end
       status = fprintf(fid,'\n%10s',' ');
       Zerr = allData{j}.Zerr;
       if allData{j}.Cmplx
           for i = 1:nComp/2
               status = fprintf(fid,'%15.6E %15.6E',Zerr(k,i),Zerr(k,i));
           end
       else
           for i = 1:nComp
               status = fprintf(fid,'%15.6E',Zerr(k,i));
           end
       end
       status = fprintf(fid,'\n');
   end

end

% for plotting convenience, write origin lat & lon, if these exist
if isfield(allData{1},'origin')
    status = fprintf(fid,'\nOrigin: %12.6f %12.6f %12.6f\n',allData{1}.origin);
end

fclose(fid);
