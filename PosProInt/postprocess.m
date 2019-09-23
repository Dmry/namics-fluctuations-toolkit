%% Simulation parameters
nx=50;
ny=50;
nz=30;
dx=1;
component=0;
filename = 'mesodyn-1561646813_';
nfiles=5;

%% Function Calls
hqa=zeros(nx,1);
for filenumber=1:nfiles
    
    fp=sprintf('%s%d.vtk',filename,filenumber);
    [rho]=readvtk(component,nx,ny,nz,dx,fp);
    [rgibbs,hx,hy,q,hqx,hqy]=calculateparameters(rho,nx,ny,nz,dx);
    hq=sqrt(hqx.^2 +hqy.^2);
    hqa=hq+hqa;
    
end

hq=hqa/nfiles;
figure(1)
loglog(q,hq,'-r*');
xlabel('q')
ylabel('|h(q)|^2 A');