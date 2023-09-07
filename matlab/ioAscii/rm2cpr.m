[x,y,z,rho,nzAir] = read_mackie3d_model('Test.rm')
Config3D
m0 = mkCond3D(conf)
m0.v = log(1./rho);
m0.grid.dx = x;
m0.grid.dy = y;
m0.grid.dz = z;
m0.grid.Nx = length(x);
m0.grid.Ny = length(y);
m0.grid.NzEarth = length(z);
m0.grid.NzAir = nzAir;
[status] = writeCond_3D('TestForGary.cpr',m0,0)
