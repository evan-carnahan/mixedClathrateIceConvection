function [] = plot_flownet(Nh,Ns,h,PSI,h_style,psi_style,Grid)
   % author: Evan Carnahan
   % date: 10/20/2019
   % Description: Plots a flownet with an Nh equally spaced head contours
   %              and Ns equally spaced streamlines.
   % Input: Nh = number of head contours
   %        Ns = number of streamlines
   %        h = Grid.N by 1 column vector of heads
   %        PSI = Ny+1 by Nx+1 matrix containing the stream function
   %        h_style   = string specifying the linestyle for the head contours
   %        psi_style = string specifying the linestyle for the streamlines
   %        Grid = structure containing all information about the grid.
   % Output: none
    cStrVal = linspace(min(PSI(:)),max(PSI(:)),Ns);
    cHVal = linspace(min(h(:)),max(h(:)),Nh);
    %% plot numerical and analytical solution
    figure; 
    set(gcf,'Position',[100 100 800 200]);
    hold on
    contour(Grid.xc,Grid.yc,reshape(h,[Grid.Ny,Grid.Nx]),h_style,'LevelList',cHVal);
    contour(Grid.xf,Grid.yf,PSI,psi_style,'LevelList',cStrVal);
    %daspect([1 1 1]);
    xlabel('x-dir, m')
    ylabel('z-dir, m')