i=0;
inversion = 'Modular_NLCG';
loc = 'default/';
modelFile = [loc inversion '_' sprintf('%03d',i) '.rho'];
dataFile = [loc inversion '_' sprintf('%03d',i) '.dat'];
[AllData] = readZ_3D(dataFile,'[mV/km]/[nT]');
while exist(modelFile,'file')
    [Cond] = readCond_3D(modelFile,2);
    if strcmp(Cond.paramType,'LOGE')
        Cond.paramType = 'LOG10';
        Cond.v = Cond.v / log(10);
        Cond.AirCond = Cond.AirCond / log(10);
    elseif strcmp(Cond.paramType,'LINEAR')
        Cond.paramType = 'LOG10';
        Cond.v = log10(Cond.v);
        Cond.AirCond = log10(Cond.AirCond);
    end
    
    plot3Dview
    caxis([-3.5 -0.5]);
    %cb=colorbar('ytick',[-3,-2,-1],'yticklabel',{'1000','100','10'},'FontSize',...
    %    16,'FontWeight','demi');
    %title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
    %shading flat
    
    t=title(['NLCG iteration: ' num2str(i)],'FontWeight','demi','FontSize',24);
    %set(t,'position',[-115 45 -80]);
    % title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
    % set(cb,'ytick',[-3 -2.5 -2 -1.5 -1 -0.5 0])
    % set(cb,'yticklabel',[1000 300 100 30 10 3 1],'FontWeight','demi','FontSize',18)
    % set(cb,'position',get(cb,'position')+[0.05 0 0 0])
    % view([-40,15])
    print('-dpng','-r300',[loc inversion '_' sprintf('%03d',i) '.png']);
    
    fclose all; close all;
    i = i+1;
    modelFile = [loc inversion '_' sprintf('%03d',i) '.rho'];
end

