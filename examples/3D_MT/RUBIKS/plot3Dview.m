%   sliced open 3D plots for ModEM paper:
%    using condPlot to load a conductivity file: everything needed is in
%     data structure Cond + load data
%%
x = cumsum([0; Cond.grid.dx])+Cond.grid.origin(1);
y = cumsum([0; Cond.grid.dy])+Cond.grid.origin(2);
z = cumsum([0; Cond.grid.dz]);
[Xxy Yxy] = meshgrid(x,y);
[Xxz Zxz] = meshgrid(x,z);
[Yyz Zyz] = meshgrid(y,z);
%   pad cond
N = size(Cond.v);
v = zeros(N+1);
v(1:end-1,1:end-1,1:end-1) = Cond.v;

%   block 1: 
%   surface 1: part  of horizontal plane (level iz = 5 ... 2km)
%    iy 2:end/2-1; iy = 3:end-2
figure('position',[100,100,800,600],'PaperPosition',[1,1,8,6]);
iz = 5; 
ix1 = 3; ix2 = N(1)-1;
iy1 = 3; iy2 = round(N(2)/2);
surf(Xxy(iy1:iy2,ix1:ix2)',Yxy(iy1:iy2,ix1:ix2)', ...
    z(iz)*ones(size(Xxy(ix1:ix2,iy1:iy2))),v(ix1:ix2,iy1:iy2,iz))
set(gca,'Zdir','reverse','Xdir','reverse','Ydir','reverse')
caxis([-3.5,-.5])
hold on
%   surface 2: yz plane: same range of y, z from iz to  layer 45
iz1 = iz; iz2 = 41;
surf(x(ix2)*ones(size(Yyz(iz1:iz2,iy1:iy2)))',Yyz(iz1:iz2,iy1:iy2)', ...
    Zyz(iz1:iz2,iy1:iy2)',squeeze(v(ix2,iy1:iy2,iz1:iz2)));
%   surface 3: xz plane: same range of x, z from iz to  layer 45
surf(Xxz(iz1:iz2,ix1:ix2)',y(iy2)*ones(size(Xxz(iz1:iz2,ix1:ix2)))', ...
    Zxz(iz1:iz2,ix1:ix2)',squeeze(v(ix1:ix2,iy2,iz1:iz2)));

%   block 2: 
iz = 5; 
ix1 = 3; ix2 = round(N(1)/2);
iy1 = iy2; iy2 = round(N(2)/2)+18;
surf(Xxy(iy1:iy2,ix1:ix2)',Yxy(iy1:iy2,ix1:ix2)', ...
    z(iz)*ones(size(Xxy(ix1:ix2,iy1:iy2))),v(ix1:ix2,iy1:iy2,iz))
set(gca,'Zdir','reverse','Xdir','reverse','Ydir','reverse')
caxis([-3.5,-.5])
hold on
%   surface 2: yz plane: same range of y, z from iz to  layer 45
surf(x(ix2)*ones(size(Yyz(iz1:iz2,iy1:iy2)))',Yyz(iz1:iz2,iy1:iy2)', ...
    Zyz(iz1:iz2,iy1:iy2)',squeeze(v(ix2,iy1:iy2,iz1:iz2)));
%   surface 3: xz plane: same range of x, z from iz to  layer 45
surf(Xxz(iz1:iz2,ix1:ix2)',y(iy2)*ones(size(Xxz(iz1:iz2,ix1:ix2)))', ...
    Zxz(iz1:iz2,ix1:ix2)',squeeze(v(ix1:ix2,iy2,iz1:iz2)));

%stcor = AllData{1}.siteLoc/1000;
%plot3(stcor(:,1),stcor(:,2),stcor(:,3),'ko','MarkerSize',7,'MarkerFaceColor',...
%    [1,1,1],'Linewidth',2);
set(gca,'Fontsize',16,'FontWeight','demi')
zlabel('km')
%cb=colorbar('ytick',[-3,-2,-1],'yticklabel',{'1000','100','10'},'FontSize',...
%    16,'FontWeight','demi');
%title(cb,'\rho [{\Omega}m]','FontWeight','demi','FontSize',18)
%posn = get(cb,'position');
%set(cb,'position',[posn(1) 3*posn(2) 2*posn(3) posn(4)/2]);
%shading flat