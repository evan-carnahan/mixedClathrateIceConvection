% wrapper_stokes.m
% authors: Marc Hesse, Evan Carnahan
clc, close all, clear
% Outline of pressure grid
Gridp.xmin = 0; Gridp.xmax = 1; Gridp.Nx = 100;
Gridp.ymin = 0; Gridp.ymax = 1; Gridp.Ny = 100;

Grid = build_stokes_grid(Gridp);
[D,Edot,Dp,Gp,I]=build_stokes_ops(Grid);
A = D*Edot;
%incompressible
Zp = zeros(Grid.p.N);
%compressible

% 
% subplot 131
% spy(Edot)
% subplot 132
% spy(D)
% subplot 133
% spy(A-A')

L = [A,Gp;
    Dp,Zp];

%need to elminate corners
Param.dof_dir =  [...
                  Grid.x.dof_ymax;...           %set y_max x-vel
                  Grid.x.dof_xmax;...           %set x_max x-vel
                  Grid.x.dof_xmin;...           %set x_min x-vel
                  Grid.x.dof_ymin;...           %set y_min x-vel
                  Grid.x.N+Grid.y.dof_ymin;...  %set y_min y-vel
                  Grid.x.N+Grid.y.dof_ymax;...  %set y_max y-vel
                  Grid.x.N+Grid.y.dof_xmin;...  %set x_min y-vel
                  Grid.x.N+Grid.y.dof_xmax;...  %set x_max y-vel
                  Grid.p.Nf+Grid.p.Ny+Grid.p.Nx+2%set pressure
                  ];                            
Param.g = 0*Param.dof_dir;
Param.g(1:(Grid.x.Nx),1) = 1;
%Param.dor_dir_uni = unique(Param.dof_dir)
B = I([Param.dof_dir],:);
N = I;
N(:,[Param.dof_dir]) = [];
fs = 0;
Lr = N'*L*N;
find(all(full(Lr)==0,2))
u = solve_lbvp(L,fs,B,Param.g,N);

ux = reshape(u(1:Grid.x.N),Grid.x.Ny,Grid.x.Nx)
uy = reshape(u(Grid.p.Nfx+1:Grid.p.Nf),Grid.p.Ny+1,Grid.p.Nx)
up = reshape(u(Grid.p.Nf+1:end),Grid.p.Ny,Grid.p.Nx)


[ux_x,ux_y] = meshgrid(Grid.p.xf,Grid.p.yc);
[uy_x,uy_y] = meshgrid(Grid.p.xc,Grid.p.yf);
% uy_y = flip(uy_y);
% ux_y = flip(ux_y);
[up_x,up_y] = meshgrid(Grid.p.xc,Grid.p.yc);
% up_y = flip(up_y);
ux_pVal = interp2(ux_x,ux_y,ux,up_x,up_y);
uy_pVal = interp2(uy_x,uy_y,uy,up_x,up_y);

figure;
hold on
quiver(ux_x,ux_y,ux,zeros(size(ux)),'b')
quiver(uy_x,uy_y,zeros(size(uy)),uy,'g')
plot(ux_x,ux_y,'ob')
plot(uy_x,uy_y,'og')
title('x-y velocity on calculated grid')

figure;
hold on
quiver(up_x,up_y,ux_pVal,uy_pVal,'r')
plot(up_x,up_y,'or')
title('Interpolated velocity')

figure;
hold on
contourf(up_x,up_y,up-mean(up)); colorbar()
clim((-1,1));
% plot(up_x,up_y,'or')
title('Pressure')

figure;
q = u(1:Grid.p.Nf);
[Xp,Yp] = meshgrid(Grid.p.xf,Grid.p.yf);
Qymin = [0 cumsum(q(Grid.p.dof_f_ymin))'];      % Integral of flow int
Qx    = reshape(q(1:Grid.p.Nfx),Grid.p.Ny,Grid.p.Nx+1);  % Horizontal fluxes
PSI   = cumsum([-Qymin;Qx],1); % integrals into domain with Qymin as initial value
surf(Xp,Yp,PSI);



%% For next time: Add in pressure change-> compressible stokes
%-> Possibly talk to Jake on pressure block
%-> Look at corner fluxes velocity split in opposite direction along middle
%       this will allow for a check for analytical solution -> Email Marc
%       for details






