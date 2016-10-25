clear
close
clc;
nx=10;
ny=10;

count=0;
tfinal=1;
xlength=1;
ylength=1;
dx=xlength/(nx-1);
dy=ylength/(ny-1);

disp(' ')
disp(' Explicit Time Stepping, 2D Heat Equation')
disp(' ')

alpha=0.25;
dt1=alpha*dx^2;
dt2=alpha*dy^2;
dt=min(dt1,dt2);

for i=1:nx;
    for j=1:ny;
        T(i,j)=0;
    end
end
% start values 0.d0 on x for y=0 
for i=1:nx;
    T(i,1)=0;
end
% start values 1.d0 on x for y=1 
for i=1:nx;
    T(i,ny)=1;
end
% start values 1.d0 on y for x=0 
for j=2:ny-1;
    T(1,j)=T(1,j-1) + dy;
end
% start values 1.d0 on y for x=1 
for j=2:ny-1;
    T(nx,j)=T(nx,j-1) + dy;
end

imagesc(T);     
axis xy;   colorbar;  axis square;  axis image;
disp('pause'), disp(' '), pause

for i=1:nx;
    for j=1:ny;
        Tnew(i,j)=T(i,j);
    end
end

dphimax=1.0e+08; tolerance=1.0e-05;
time=0; itmax = 15000; it = 0;
% explicit time iteration 
tstart = tic;
while (dphimax>=tolerance & it < itmax);
    it = it+1
    time=time +dt;
    dphimax = 0; 
    for i=2:nx-1;
        for j=2:ny-1;
            dphi=dt*(T(i+1,j)-2*T(i,j)+T(i-1,j))/dx^2+ ...
                      dt*(T(i,j+1)-2*T(i,j)+T(i,j-1))/dy^2;
            Tnew(i,j)=T(i,j) + dphi;
            % max update
            dphimax = max(dphimax, dphi); 
        end
    end
    for i=1:nx;
        for j=1:ny;
            T(i,j)=Tnew(i,j);
        end
    end 
%    it
%    imagesc(T);      
%    axis xy;   colorbar;  axis square;  axis image;
%    disp('pause'), disp(' '), pause
end
telapsed = toc(tstart);
fprintf('[Explizit Time Stepping] in %d-iteration\n',it); 
fprintf('[Explizit Time Stepping] in %f sec.\n', telapsed );
    figure; 
imagesc(T); 
axis xy;   colorbar;  axis square;  axis image;
