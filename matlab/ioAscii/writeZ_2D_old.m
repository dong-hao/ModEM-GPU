function [status] = writeZ_2D(cfile,allData,info,units,isign)
%  Usage:  [status] = writeZ_2D(cfile,allData,info,units,isign);
%   arguments info, units and isign are optional.
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   NB  Assuming complex data
%  (c) Anna Kelbert, 2008

nTx = length(allData);
fid = fopen(cfile,'w');

%  description <= 80 char in length
if nargin < 3
    info = 'Synthetic 2D MT data written in Matlab';
end

%  assume SI units by default; alternative is [mV/km]/[nT] 
if nargin < 4
    units = '[V/m]/[T]';
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

   %  record # 1 in transmitter block: period, mode, number of sites
   T = allData{j}.T;
   MODE = allData{j}.Mode;
   nSites = size(allData{j}.siteLoc);
   nComp = 2;
   status = fprintf(fid,'%12.6E %5s %5d\n',T,MODE,nSites(1));

   %  record # 2 in transmitter block: site locations
   siteLoc = allData{j}.siteLoc;
   for i = 1:2
       for k = 1:nSites
           status = fprintf(fid,'%12.3f',siteLoc(k,i));
       end
       status = fprintf(fid,'\n');
   end
    
   %  header of the data & error block: will in future be more informative
   status = fprintf(fid,'%12s Re %12s Im\n',' ',' ');
   
   %  records # 3 & #4 in transmitter block: data and error bars
   for k = 1:nSites
       if isfield(allData{j},'siteChar')
           siteChar = allData{j}.siteChar(k,:);
       else
           siteChar = num2str(k);
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
fclose(fid);
