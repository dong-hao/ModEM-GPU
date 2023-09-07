%This script plots the contour maps of observed and predicted responses for
%each frequency. The user need to input the original data file and the
%response file.

clear all

fresp_name = 'caswithprior_resp.02_01'
orgdata_fname = 'casnewcorrect.data'  
rname = 'cascad'

%Reading the response file

fid=fopen(fresp_name,'r');
ns=fscanf(fid,'%d',1);
nf=fscanf(fid,'%d',1);
nresp=fscanf(fid,'%d\n',1);
tline=fgetl(fid);
sx=fscanf(fid,'%f\n',ns); % read in N/S site locations
tline=fgetl(fid);
sy=fscanf(fid,'%f',ns); % read in E/W site locations
f=[]; Zp= zeros(nf,nresp,ns)
for i=1:nf
    f=[f fscanf(fid,'%*s%f',1)];
    temp=fscanf(fid,'%f',nresp*ns);
    Zp(i,:,:)=reshape(temp,nresp,ns);
end
fclose(fid);

%%%% Added to compute Rho and Phase  %%%%

k=796;

Zxxp = complex(Zp(:,1,:),Zp(:,2,:)).*k;
Zxyp = complex(Zp(:,3,:),Zp(:,4,:)).*k;
Zyxp = complex(Zp(:,5,:),Zp(:,6,:)).*k;
Zyyp = complex(Zp(:,7,:),Zp(:,8,:)).*k;

for i=1:nf
    temp1 = 0.2 * f(i).* (abs(Zxxp(i,:,:))).^2;
    Rhoxxp(i,:,:) = temp1;
    temp2 = 0.2 * f(i).* (abs(Zxyp(i,:,:))).^2;
    Rhoxyp(i,:,:) = temp2;
    temp3 = 0.2 * f(i).* (abs(Zyxp(i,:,:))).^2;
    Rhoyxp(i,:,:) = temp3;
    temp4 = 0.2 * f(i).* (abs(Zyyp(i,:,:))).^2;
    Rhoyyp(i,:,:) = temp4;

    Phsxxp(i,:,:) = 180.0/pi * atan2(-imag(Zxxp(i,:,:)),real(Zxxp(i,:,:)));
    Phsxyp(i,:,:) = 180.0/pi * atan2(-imag(Zxyp(i,:,:)),real(Zxyp(i,:,:)));
    Phsyxp(i,:,:) = 180.0/pi * atan2(-imag(Zyxp(i,:,:)),real(Zyxp(i,:,:)));
    Phsyyp(i,:,:) = 180.0/pi * atan2(-imag(Zyyp(i,:,:)),real(Zyyp(i,:,:)));
end

sx=sx./1000; sy=sy./1000;

% Now we read in original data responses, imposed errors, and error map

fid=fopen(orgdata_fname,'r');
tline=fgetl(fid)
while (isempty(strfind(tline,'DATA'))) % Find the starting point
    tline=fgetl(fid);
end
Zm=[];
temp=fscanf(fid,'%f\n',nresp*ns);
Zm(1,:,:)=reshape(temp,nresp,ns);
for j=2:nf
    tline=fgetl(fid);
    temp=fscanf(fid,'%f\n',nresp*ns);
    Zm(j,:,:)=reshape(temp,nresp,ns);
end
Err=[]; % reading in imposed errors
for j=1:nf
    fgetl(fid);
    temp=fscanf(fid,'%f\n',nresp*ns);
    Err(j,:,:)=reshape(temp,nresp,ns);
end
Errmap=[]; % reading in the error map
for j=1:nf
    fgetl(fid);
    temp=fscanf(fid,'%f\n',nresp*ns);
    Errmap(j,:,:)=reshape(temp,nresp,ns);
end
fclose(fid);

Zxxo = complex(Zm(:,1,:),Zm(:,2,:)).*k;
Zxyo = complex(Zm(:,3,:),Zm(:,4,:)).*k;
Zyxo = complex(Zm(:,5,:),Zm(:,6,:)).*k;
Zyyo = complex(Zm(:,7,:),Zm(:,8,:)).*k;

for i=1:nf
    temp1 = 0.2 * f(i).* (abs(Zxxo(i,:,:))).^2;
    Rhoxxo(i,:,:) = temp1;
    temp2 = 0.2 * f(i).* (abs(Zxyo(i,:,:))).^2;
    Rhoxyo(i,:,:) = temp2;
    temp3 = 0.2 * f(i).* (abs(Zyxo(i,:,:))).^2;
    Rhoyxo(i,:,:) = temp3;
    temp4 = 0.2 * f(i).* (abs(Zyyo(i,:,:))).^2;
    Rhoyyo(i,:,:) = temp4;

    Phsxxo(i,:,:) = 180.0/pi * atan2(-imag(Zxxo(i,:,:)),real(Zxxo(i,:,:)));
    Phsxyo(i,:,:) = 180.0/pi * atan2(-imag(Zxyo(i,:,:)),real(Zxyo(i,:,:)));
    Phsyxo(i,:,:) = 180.0/pi * atan2(-imag(Zyxo(i,:,:)),real(Zyxo(i,:,:)));
    Phsyyo(i,:,:) = 180.0/pi * atan2(-imag(Zyyo(i,:,:)),real(Zyyo(i,:,:)));
end

close all;
% Plotting starts here
% Plotting Phs YX predicted component

for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsyxp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 -90])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs YX predicted plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['PhsYXp',num2str(round(f(i))),'s'])
end
% Plotting Phs YX observed component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsyxo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 -90])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs YX observed plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['PhsYXo',num2str(round(f(i))),'s'])
end
% Plotting Phs XY observed component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsxyo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 90])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs XY observed plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['PhsXYo',num2str(round(f(i))),'s'])
end
% Plotting Phs XY predicted component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsxyp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 90])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs XY predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['PhsXYp',num2str(round(f(i))),'s'])
end
% Plotting Phs XX observed
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsxxo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 180])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs XX observed plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['PhsXXobserved_',num2str(round(f(i))),'s'])
end
% Plotting Phs XX predicted
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsxxp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 180])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs XX predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['PhsXXpredicted_',num2str(round(f(i))),'s'])
end
% Plotting Phs YY observed
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsyyo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 180])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs YY observed plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['PhsYYobserved_',num2str(round(f(i))),'s'])
end
% Plotting Phs YY predicted
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Phsyyp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-180 180])
    xlabel(xlab)
    ylabel(ylab)
    title(['Phs YY predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['PhsYYpredicted_',num2str(round(f(i))),'s'])
end
% Plotting Rho YX predicted component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = log10(Rhoyxp(i,:));
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 2.5])
    xlabel(xlab)
    ylabel(ylab)
    title(['Rho YX predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['RhoYXp',num2str(round(f(i))),'s'])
end
% Plotting Rho YX observed component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = log10(Rhoyxo(i,:));
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 2.5])
    xlabel(xlab)
    ylabel(ylab)
    title(['Rho YX observed plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['RhoYXo',num2str(round(f(i))),'s'])
end
% Plotting Rho XY observed component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = log10(Rhoxyo(i,:));
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 2.5])
    xlabel(xlab)
    ylabel(ylab)
    title(['Rho XY observed plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['RhoXYo',num2str(round(f(i))),'s'])
end
% Plotting Rho XY predicted component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = log10(Rhoxyp(i,:));
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[0 2.5])
    xlabel(xlab)
    ylabel(ylab)
    title(['Rho XY predicted plot for ',' ',num2str(f(i)),' s'])
    colorbar
    print('-depsc2',['RhoXYp',num2str(round(f(i))),'s'])
end
% Plotting Rho XX observed component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Rhoxxo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-20 20])
    xlabel(xlab)
    ylabel(ylab)
    title(['RhoXX observed plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['RhoXXobserved_',num2str(round(f(i))),'s'])
end
% Plotting Rho XX predicrted component
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Rhoxxp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-20 20])
    xlabel(xlab)
    ylabel(ylab)
    title(['RhoXX predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['RhoXXpredicted_',num2str(round(f(i))),'s'])
end
% Plotting Rho YY observed
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Rhoyyo(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-20 20])
    xlabel(xlab)
    ylabel(ylab)
    title(['RhoYY observed plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['RhoYYobserved_',num2str(round(f(i))),'s'])
end
% Plotting Rho YY predicted
for i=1:nf
    figure (i)
    colormap(flipud(jet(32)));
    x = (-420:10:450);
    xlims=[-420 450]; ylims=[-420 420];
    xlab='distance in EW direction (km)';
    ylab='distance in NS direction (km)';
    value = Rhoyyp(i,:);
    [x,y]=meshgrid(x,x);
    [x,y,v] = griddata(sx,sy,value,x,y);
    pcolor(y,x,v);
    shading flat
    hold on
    plot(sy,sx,'k^')
    axis equal
    set(gca,'xlim',xlims,'ylim',ylims)
    set(gca,'clim',[-20 20])
    xlabel(xlab)
    ylabel(ylab)
    title(['RhoYY predicted plot for ',' ',num2str(f(i)),' s'])

    colorbar
    print('-depsc2',['RhoYYpredicted_',num2str(round(f(i))),'s'])
end

