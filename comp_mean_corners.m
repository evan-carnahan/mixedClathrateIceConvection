function [Kd] = comp_mean_corners(K,p,Grid)
   % author: Evan Carnahan
   % date: 10/21/2019
   % Description:
   % Takes coefficient field, K, defined at the cell centers and computes the
   % mean specified by the power, p and returns it in a sparse diagonal
   % matrix, Kd.
   %
   % Input:
   % K = Ny by Nx matrix of cell centered values
   % p = power of the generalized mean
   %       1 (arithmetic mean)
   %      -1 (harmonic mean)
   % Grid = structure containing information about the grid.
   %
   % Output:
   % Kd = Nf by Nf diagonal matrix of power means at the cell faces.
   %
   % Example call:
   % K = @(x) 1+x.^3;
   % Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
   % Grid = build_grid(Grid);
   % Kd = comp_mean(K(Grid.xc),1,Grid);
   
% mean_k = zeros(Grid.Ny+1,Grid.Nx+1);
% mean_k(2:Grid.Ny,2:Grid.Nx) = 4*(K(1:Grid.Ny-1,1:Grid.Nx-1).^p+K(2:Grid.Ny,2:Grid.Nx).^p...
%     +K(1:Grid.Ny-1,2:Grid.Nx).^p+K(2:Grid.Ny,1:Grid.Nx-1).^p).^(1/p);
% 
% mean_k(Grid.Ny+1,2:Grid.Nx) = 2*(K(Grid.Ny,1:Grid.Nx-1).^p+K(Grid.Ny,2:Grid.Nx).^p).^(1/p);
% mean_k(1,2:Grid.Nx) = 2*[K(1,2:Grid.Nx).^p+K(1,1:Grid.Nx-1).^p].^(1/p);
% 
% mean_k(2:Grid.Ny,Grid.Nx+1) = 2*(K(2:Grid.Ny,Grid.Nx).^p+K(1:Grid.Ny-1,Grid.Nx).^p).^(1/p);
% mean_k(2:Grid.Ny,1) = 2*(K(1:Grid.Ny-1,1).^p+K(2:Grid.Ny,1).^p).^(1/p);

mean_k = nan(Grid.Ny+1,Grid.Nx+1);
% mean_k = zeros(Grid.Ny+1,Grid.Nx+1);
mean_k(2:Grid.Ny,2:Grid.Nx) = ((K(1:Grid.Ny-1,1:Grid.Nx-1).^p+K(2:Grid.Ny,2:Grid.Nx).^p...
    +K(1:Grid.Ny-1,2:Grid.Nx).^p+K(2:Grid.Ny,1:Grid.Nx-1).^p)/4).^(1/p);

mean_k(Grid.Ny+1,2:Grid.Nx) = ((K(Grid.Ny,1:Grid.Nx-1).^p+K(Grid.Ny,2:Grid.Nx).^p)/2).^(1/p);

mean_k(1,2:Grid.Nx) = ((K(1,2:Grid.Nx).^p+K(1,1:Grid.Nx-1).^p)/2).^(1/p);

mean_k(2:Grid.Ny,Grid.Nx+1) = ((K(2:Grid.Ny,Grid.Nx).^p+K(1:Grid.Ny-1,Grid.Nx).^p)/2).^(1/p);
mean_k(2:Grid.Ny,1) = ((K(1:Grid.Ny-1,1).^p+K(2:Grid.Ny,1).^p)/2).^(1/p);

KdVec = mean_k(:);
Kd = KdVec;
%Kd = spdiags(mean_k(:),0,(Grid.Nx+1)*(Grid.Ny+1),(Grid.Nx+1)*(Grid.Ny+1));

   