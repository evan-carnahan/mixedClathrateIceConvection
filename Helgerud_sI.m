function rho = Helgerud_sI(P,T)%out = Helgerud_sI(P,T)
    % from Helgerud et al. (2009) - I think for 97.6% cage filling
    % P in MPa
    % T in K then converted to C
    T=T-273.15;
%     out.Vp = -1.84.*T + 0.31.*P + 3766;
%     out.Vs = -0.892.*T -0.1.*P + 1957;
%     out.shear = -4.2e-3.*T + 9e-5.*P + 3.541;
%     out.K = -1.9e-2.*T + 3.8e-3.*P + 8.39;
%     out.rho = (-2.3815e-4.*T + 1.1843e-4.*P +0.91801)*1e3;
rho = (-2.3815e-4.*T + 1.1843e-4.*P +0.91801)*1e3;