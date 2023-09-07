function [allData,header,units,isign,origin,info] = readZ_3D(cfile,newunits,onetype)
%  Usage:  [allData,header,units,isign,origin,info] = readZ_3D(cfile,newunits,onetype);
%   read contents of cell array allData from file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%   If onetype is specified, only read that data type and skip others.
%  (c) Anna Kelbert, 2011-2013


fid = fopen(cfile,'r');

if nargin < 3
    onetype = '';
end

%  read the data: one block per data type
while 1
    % read the header info
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    tmp = char(tmp{1});
    if isempty(tmp); break; end
    header = tmp;
    % read the block header
    tmp = textscan(fid,'# %s',1,'delimiter','\n');
    blockinfo = textscan(fid,'> %s',6,'delimiter','\n');
    blockinfo = char(blockinfo{1});
    dataType = blockinfo(1,:);
    signstr  = blockinfo(2,:);
    if findstr(signstr,'-')
        isign = -1;
    else
        isign = +1;
    end
    typeUnits = strtrim(blockinfo(3,:));
    orientation = sscanf(blockinfo(4,:),'%f',1);
    origin = sscanf(blockinfo(5,:),'%f',2);
    tmp  = sscanf(blockinfo(6,:),'%d',2);
    nTx = tmp(1);
    nSites = tmp(2);
    % now read the data line by line
    switch strtrim(dataType)
        case 'Full_Impedance'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                continue
            end
            if nargin > 1
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 4;
            comp = ['ZXX';'ZXY';'ZYX';'ZYY'];
            info{1}.data = nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{1}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{1}.err(k,j,i) = data{11}(ind);
                            info{1}.lat(k) = data{3}(ind);
                            info{1}.lon(k) = data{4}(ind);
                            info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.ncomp = 8;
            info{1}.comp = comp;
        case 'Off_Diagonal_Impedance'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                continue
            end
            if nargin > 1
                SI_factor = ImpUnits(typeUnits,newunits);
                units = newunits;
            else
                SI_factor = 1.0;
                units = typeUnits; % output the units of impedances
            end
            ncomp = 2;
            comp = ['ZXY';'ZYX'];
            info{1}.data = nan(nSites,nTx,ncomp);
            info{1}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(strtrim(char(unique(data{2}))));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{1}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{1}.err(k,j,i) = data{11}(ind);
                            info{1}.lat(k) = data{3}(ind);
                            info{1}.lon(k) = data{4}(ind);
                            info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.ncomp = 4;
            info{1}.comp = comp;
        case 'Full_Vertical_Components'
            if ~isempty(onetype) && ~strcmp(onetype,strtrim(dataType))
                continue
            end
            ncomp = 2;
            comp = ['TX ';'TY '];
            info{2}.data = nan(nSites,nTx,ncomp);
            info{2}.err = nan(nSites,nTx,ncomp);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                codelist = strtrim(char(data{2}(itx)));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(strtrim(codes(k,:)),codelist);
                    complist = strtrim(char(data{8}(itx(irx))));
                    for i = 1:ncomp % ... and store all components
                        icomp = strmatch(strtrim(comp(i,:)),complist);
                        if ~isempty(icomp)
                            ind = itx(irx(icomp));
                            info{2}.data(k,j,i) = data{9}(ind) + 1i* data{10}(ind);
                            info{2}.err(k,j,i) = data{11}(ind);
                            info{2}.lat(k) = data{3}(ind);
                            info{2}.lon(k) = data{4}(ind);
                            info{2}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                        end
                    end
                end
            end
            info{2}.code = codes;
            info{2}.per = periods;
            info{2}.ncomp = 4;
            info{2}.comp = comp;
        otherwise
            disp('Unknown data type');
            break;
    end   
end
status = fclose(fid);

% if two data types are present, merge the periods
if length(info) == 1
    per = info{1}.per;
    ind1 = 1:length(per);
else
    per = sort(union(info{1}.per,info{2}.per));
    [tmp,ind1] = intersect(per,info{1}.per);
    [tmp,ind2] = intersect(per,info{2}.per);
end

% compute the maximum total number of components
ncomp = 0;
for i = 1:length(info)
    ncomp = ncomp + info{i}.ncomp/2;
end

% for compatibility, convert to the (old) allData structure - most programs
% expect full impedances and vertical components (in this order) 
for j = 1:length(per)
    for k = 1:length(info)
        ind = find(info{k}.per==per(j));
        if ~isempty(ind); itx(k) = ind; else itx(k) = 0; end
    end    
    allsites = info{1}.code;
    allsitesloc = info{1}.loc;
    allsiteslat = info{1}.lat;
    allsiteslon = info{1}.lon;
    for k = 2:length(info)
        for i = 1:length(info{k}.code)
            if isempty(intersect(allsites,info{k}.code(i,:),'rows'))
                allsites = [allsites; info{k}.code(i,:)];
                allsitesloc = [allsitesloc; info{k}.loc(i,:)];
                allsiteslat = [allsiteslat; info{k}.lat(i)];
                allsiteslon = [allsiteslon; info{k}.lon(i)];
            end
        end
    end
    nsites = length(allsites);
    allData{j}.T = per(j);
    allData{j}.Cmplx = 1;
    allData{j}.units = units;
    allData{j}.signConvention = isign;
    allData{j}.nComp = ncomp;
    %allData{j}.compChar(1:ncomp) = '';
    allData{j}.siteLoc = allsitesloc;
    allData{j}.siteChar = allsites;
    allData{j}.Z(1:nsites,1:ncomp) = NaN;
    allData{j}.Zerr(1:nsites,1:ncomp) = NaN;
    allData{j}.nComp = 2*ncomp;
    allData{j}.origin = [origin' 0];
    allData{j}.orient = orientation;
    allData{j}.lat = allsiteslat';
    allData{j}.lon = allsiteslon';
    icomp1 = 1;
    for k = 1:length(info)
        if itx(k) > 0
            icomp2 = icomp1 + info{k}.ncomp/2 - 1;
            [sites,irx] = intersect(allsites,info{k}.code,'rows');
            allData{j}.Z(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.data(:,itx(k),:));
            allData{j}.Zerr(irx,icomp1:icomp2) = SI_factor*squeeze(info{k}.err(:,itx(k),:));
            allData{j}.compChar(icomp1:icomp2,:) = info{k}.comp;
            icomp1 = icomp2+1;
        end
    end
end