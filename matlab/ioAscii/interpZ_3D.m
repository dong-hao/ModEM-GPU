function newdata = interpZ_3D(data,siteLoc,siteChar,pererr)

% newdata = interpZ_3D(data,siteLoc,siteChar,pererr)
%
% accepts a 3D MT data structure and interpolates to user-specified
% data locations; saves in a new data structure
%
% optionally replaces error bars with pererr percent (e.g. pererr = 5)
%
% use in conjuction with
% data = readZ_3D(dataFile);
% writeZ_3D(dataFileErrors,data);

if nargin > 3
    ADD_NOISE = 1;
else
    ADD_NOISE = 0;
end

newdata = data;
newx = siteLoc(:,1);
newy = siteLoc(:,2);

for k = 1:length(data)
    nComp = size(data{k}.Z,2);
    oldx = data{k}.siteLoc(:,1);
    oldy = data{k}.siteLoc(:,2);
    olddat = data{k}.Z;
    olderr = data{k}.Zerr;
    newdat = zeros(length(newx),nComp);
    newerr = zeros(length(newx),nComp);
    for j = 1:nComp
        newdat(:,j) = griddata(oldx,oldy,olddat(:,j),newx,newy,'cubic');
        newerr(:,j) = griddata(oldx,oldy,olderr(:,j),newx,newy,'cubic');
    end    
    newdata{k}.siteLoc = siteLoc;
    newdata{k}.siteChar = siteChar;
    newdata{k}.Z = newdat;
    newdata{k}.Zerr = newerr;
end

if ADD_NOISE
    disp(['adding noise ... % Error = ' num2str(pererr)]);
    FracError = pererr/100;
    newdata = addNoise_3D(newdata,FracError);
end


