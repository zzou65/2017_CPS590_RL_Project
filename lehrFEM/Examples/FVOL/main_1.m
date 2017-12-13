% runs the exp_1.m function for solving a convection diffusion equation

% -eps grad div u + [1 0] grad u
% exact solution given (see documenatation)
% dirichlet boundary
% computational domain: triangle [0 0],[1 1],[1 -1]

% define constants

NREFS = 2;            % number of refinement steps
eps = 10^-6;          % diffusion constant
plot = 0;             % plotting flag
output = 0;           % output flag

% solve problem

[errs, mw, msh, u] =exp_1(eps,NREFS,plot,output);

% plot the solution

plot_Mesh(msh,'as')
plot_LFEfancy(u,msh);
plot_LFE(u,msh); colorbar;

disp(['MeshWidth = ',num2str(mw)]);

% clear memory

clear all;