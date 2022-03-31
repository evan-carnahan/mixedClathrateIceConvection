function [VXc, VYc] = robin_hood(vm,Grid)

[X,Y]     = meshgrid(Grid.xc,Grid.yc);
VX = reshape(vm(1:Grid.Nfx),Grid.Ny,Grid.Nx+1);
[Xx,Yx] = meshgrid(Grid.xf,Grid.yc);
VY = reshape(vm(Grid.Nfx+1:Grid.Nf),Grid.Ny+1,Grid.Nx);
[Xy,Yy] = meshgrid(Grid.xc,Grid.yf);

% Interpolate to cell centers
VXc = interp2(Xx,Yx,VX,X,Y);
VYc = interp2(Xy,Yy,VY,X,Y);


end