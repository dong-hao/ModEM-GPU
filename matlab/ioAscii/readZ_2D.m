function [allData,header,units,isign,origin,info] = readZ_2D(cfile)
%  Usage:  [allData,header,units,isign,origin,info] = readZ_2D(cfile);
%   read contents of cell array allData from file
%   cfile.  There is one cell per period; each
%   cell contains all information necessary to define
%   data (locations, values, error standard dev) for
%   each period (transmitter)
%   Site locations units have to match the model, i.e. use meters.
%  (c) Anna Kelbert, 2011

fid = fopen(cfile,'r');

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
        case 'TE_Impedance'
            info{1}.data = nan(nSites,nTx);
            info{1}.err = nan(nSites,nTx);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(codes(k,:),char(data{2}(itx)));
                    ind = itx(irx);
                    info{1}.data(k,j) = data{9}(ind) + 1i* data{10}(ind);
                    info{1}.err(k,j) = data{11}(ind);
                    info{1}.lat(k) = data{3}(ind);
                    info{1}.lon(k) = data{4}(ind);
                    info{1}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                end
            end
            info{1}.code = codes;
            info{1}.per = periods;
            info{1}.mode = 'TE';
            units = typeUnits; % output the units of impedances
        case 'TM_Impedance'
            info{2}.data = nan(nSites,nTx);
            info{2}.err = nan(nSites,nTx);
            data = textscan(fid,'%f %s %f %f %f %f %f %s %f %f %f');
            periods = unique(data{1});
            codes = sortrows(char(unique(data{2})));
            for j = 1:length(periods)
                % find all data for j'th period
                itx = find(data{1}==periods(j));
                for k = 1:length(codes)
                    % find all data for k'th site for j'th period
                    irx = strmatch(codes(k,:),char(data{2}(itx)));
                    ind = itx(irx);
                    info{2}.data(k,j) = data{9}(ind) + 1i* data{10}(ind);
                    info{2}.err(k,j) = data{11}(ind);
                    info{2}.lat(k) = data{3}(ind);
                    info{2}.lon(k) = data{4}(ind);
                    info{2}.loc(k,1:3) = [data{5}(ind) data{6}(ind) data{7}(ind)];
                end
            end
            info{2}.code = codes;
            info{2}.per = periods;
            info{2}.mode = 'TM';
        otherwise
            disp('Unknown data type');
            break;
    end
end
status = fclose(fid);

% for compatibility, convert to the (old) allData structure - assume full
% impedances and vertical components (in this order) with same periods etc
ind = 1;
for j = 1:length(info{1}.per)
    irx = find(~isnan(info{1}.data(:,j,1)));
    allData{ind}.T = info{1}.per(j);
    allData{ind}.Cmplx = 1;
    allData{ind}.Mode = info{1}.mode;
    allData{ind}.units = units;
    allData{ind}.signConvention = isign;
    allData{ind}.siteLoc = info{1}.loc(irx,2:3);
    allData{ind}.siteChar = info{1}.code(irx,:);
    allData{ind}.Z = squeeze(info{1}.data(irx,j));
    allData{ind}.Zerr = squeeze(info{1}.err(irx,j));
    allData{ind}.origin = [origin' 0];
    allData{ind}.orient = orientation;
    allData{ind}.lat = info{1}.lat(irx)';
    allData{ind}.lon = info{1}.lon(irx)';
    ind = ind+1;
end
for j = 1:length(info{2}.per)
    irx = find(~isnan(info{2}.data(:,j,1)));
    allData{ind}.T = info{2}.per(j);
    allData{ind}.Cmplx = 1;
    allData{ind}.Mode = info{2}.mode;
    allData{ind}.units = units;
    allData{ind}.signConvention = isign;
    allData{ind}.siteLoc = info{2}.loc(irx,2:3);
    allData{ind}.siteChar = info{2}.code(irx,:);
    allData{ind}.Z = squeeze(info{2}.data(irx,j));
    allData{ind}.Zerr = squeeze(info{2}.err(irx,j));
    allData{ind}.origin = [origin' 0];
    allData{ind}.orient = orientation;
    allData{ind}.lat = info{2}.lat(irx)';
    allData{ind}.lon = info{2}.lon(irx)';
    ind = ind+1;
end
