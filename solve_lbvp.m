function [u] = solve_lbvp(L,f,B,g,N)
% author: Evan Carnahan
% date: 9/20/2019
% Description
% Computes the solution $u$ to the linear differential problem given by
%
% $$\mathcal{L}(u)=f \quad x\in \Omega $$
%
% with boundary conditions
%
% $$\mathcal{B}(u)=g \quad x\in\partial\Omega$$.
%
% Input:
% L = matrix representing the discretized linear operator of size N by N,
% where N is the number of degrees of fredom
% f = column vector representing the discretized r.h.s. and contributions
% due non-homogeneous Neumann BC?s of size N by 1
% B = matrix representing the constraints arising from Dirichlet BC?s of
% sizeNcbyN
% g = column vector representing the non-homogeneous Dirichlet BC?s of size % Ncby1.
% N = matrix representing a orthonormal basis for the null-space of B and
% of size N by (N-Nc).
% Output:
% u = column vector of the solution of size N by 1
%u = N*[(N'*L*N)\N'*(f-L
u_p = B'*((B*B')\g);
u_0 = N*((N'*L*N)\(N'*(f-(L*u_p))));
u = u_0+u_p;

% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir = Grid.dof_xmin;
% >> Param.dof_f_dir = Grid.dof_f_xmin;
% >> Param.dof_neu = Grid.dof_xmax;
% >> Param.dof_f_neu = Grid.dof_f_xmax;
% >> Param.qb = 1;
% >> Param.g = 0;
% >> [B,N,fn] = build_bnd(Param,Grid,I); % Build constraint matrix and
%
% >>L=-D*G;
% >> fs = spalloc(Grid.N,1,0);
% >> h = solve_lbvp(L,fs+fn,B,Param.g,N); % Solve linear boundary value problem