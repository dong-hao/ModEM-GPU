function [status] = writeZ_2D(cfile,allData,header,units,isign)
%  Usage:  [status] = writeZ_2D(cfile,allData,header,units,isign);
%   write contents of cell array allData to file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%  (c) Anna Kelbert, 2011

%  assume 1 complex component for either TE or TM mode
TEind(1:length(allData)) = 0;
TMind(1:length(allData)) = 0;
for i = 1:length(allData)
    if strcmp(allData{i}.Mode,'TE')
        TEind(i) = 1;
    else
        TMind(i) = 1;
    end
end
allDataTE = allData(find(TEind));
allDataTM = allData(find(TMind));

% TE block

%  compute the total number of frequencies and sites
nTx = length(allDataTE);
allsites = allDataTE{1}.siteChar;
allsitesloc = allDataTE{1}.siteLoc;
if isfield(allDataTE{1},'lat')
    GG = 1;
    lat = allDataTE{1}.lat;
    lon = allDataTE{1}.lon;
else
    GG = 0;
end
for j = 2:nTx
    for i = 1:length(allDataTE{j}.siteChar)
        if isempty(strmatch(allDataTE{j}.siteChar(i,:),allsites))
            allsites = [allsites; allDataTE{j}.siteChar(i,:)];
            allsitesloc = [allsitesloc; allDataTE{j}.siteLoc(i,:)];
            if GG
                lat = [lat; allDataTE{j}.lat(i)];
                lon = [lon; allDataTE{j}.lon(i)];
            end
        end
    end
end
nSites = length(allsites);
if ~ GG %  if GG locations are not supplied, set them to zero
    lat(1:nSites) = 0;
    lon(1:nSites) = 0;
end

%  regroup into observatory bins
info{1}.data = nan(nSites,nTx);
info{1}.err = nan(nSites,nTx);
for j = 1:nTx
    for i = 1:length(allsites)
        [dummy,k]=intersect(allDataTE{j}.siteChar,allsites(i,:),'rows');
        if ~isempty(k)
            info{1}.code(i,:) = allsites(i,:);
            info{1}.loc(i,:) = [0 allsitesloc(i,:)];
            info{1}.lon(i) = lon(i);
            info{1}.lat(i) = lat(i);
            info{1}.data(i,j) = allDataTE{j}.Z(k);
            info{1}.err(i,j) = allDataTE{j}.Zerr(k);
            info{1}.per(j) = allDataTE{j}.T;
        end
    end
end
info{1}.mode = 'TE';

% TM block

%  compute the total number of frequencies and sites
nTx = length(allDataTM);
allsites = allDataTM{1}.siteChar;
allsitesloc = allDataTM{1}.siteLoc;
if isfield(allDataTM{1},'lat')
    GG = 1;
    lat = allDataTM{1}.lat;
    lon = allDataTM{1}.lon;
else
    GG = 0;
end
for j = 2:nTx
    for i = 1:length(allDataTM{j}.siteChar)
        if isempty(strmatch(allDataTM{j}.siteChar(i,:),allsites))
            allsites = [allsites; allDataTM{j}.siteChar(i,:)];
            allsitesloc = [allsitesloc; allDataTM{j}.siteLoc(i,:)];
            if GG
                lat = [lat; allDataTM{j}.lat(i)];
                lon = [lon; allDataTM{j}.lon(i)];
            end
        end
    end
end
nSites = length(allsites);
if ~ GG %  if GG locations are not supplied, set them to zero
    lat(1:nSites) = 0;
    lon(1:nSites) = 0;
end

%  regroup into observatory bins
info{2}.data = nan(nSites,nTx);
info{2}.err = nan(nSites,nTx);
for j = 1:nTx
    for i = 1:length(allsites)
        [dummy,k]=intersect(allDataTM{j}.siteChar,allsites(i,:),'rows');
        if ~isempty(k)
            info{2}.code(i,:) = allsites(i,:);
            info{2}.loc(i,:) = [0 allsitesloc(i,:)];
            info{2}.lon(i) = lon(i);
            info{2}.lat(i) = lat(i);
            info{2}.data(i,j) = allDataTM{j}.Z(k);
            info{2}.err(i,j) = allDataTM{j}.Zerr(k);
            info{2}.per(j) = allDataTM{j}.T;
        end
    end
end
info{2}.mode = 'TM';

% Write to file
fid = fopen(cfile,'w');

%  description <= 80 char in length
if nargin < 3
    header = 'Synthetic 2D MT data written in Matlab';
end

%  assume SI units by default as in ModEM; alternative is [mV/km]/[nT] 
if nargin < 4
    if isfield(allData{1},'units')
        units = allData{1}.units;
    else
        units = '[V/m]/[T]'
    end
end

%  sign convention is -1 by default
if nargin < 5
    if isfield(allData{1},'signConvention')
        isign = allData{1}.signConvention;
    else
        isign = -1;
    end
end
if isign == -1
    signstr = 'exp(-i\omega t)';
else
    signstr = 'exp(+i\omega t)';
end

%  get the origin from the first period
origin = [0 0 0];
if isfield(allData{1},'origin')
    origin = allData{1}.origin;
end

%  get orientation from the first period
orientation = 0.0;
if isfield(allData{1},'orient')
    orientation = allData{1}.orient;
end

for i = 1:length(info)
    if ~isempty(info{i}.per)
        
        nTx = size(info{i}.data,2);
        nSites = size(info{i}.data,1);
        
        %  repeat file header
        fprintf(fid,'# %s\n',header);
        
        %  data type header: comment, then six lines
        comment = 'Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error';
        fprintf(fid,'# %s\n',comment);
        fprintf(fid,'> %s\n',[info{i}.mode '_Impedance']);
        fprintf(fid,'> %s\n',signstr);
        fprintf(fid,'> %s\n',units);
        fprintf(fid,'> %.2f\n',orientation);
        fprintf(fid,'> %.3f %.3f\n',origin(1:2));
        fprintf(fid,'> %d %d\n',nTx,nSites);
        
        %  now write all impedances line by line
        for k = 1:nSites
            for j = 1:nTx
                if ~isnan(info{i}.data(k,j))
                    fprintf(fid,'%12.6E ',info{i}.per(j)); % transmitter
                    fprintf(fid,'%s %8.3f %8.3f ',info{i}.code(k,:),info{i}.lat(k),info{i}.lon(k)); % receiver
                    fprintf(fid,'%12.3f %12.3f %12.3f ',info{i}.loc(k,:)); % receiver x,y,z
                    fprintf(fid,'%s %15.6E %15.6E %15.6E\n',info{i}.mode,real(info{i}.data(k,j)),imag(info{i}.data(k,j)),info{i}.err(k,j)); % data
                end
            end
        end
        
    end
end

status = fclose(fid);
