% runs the exp_3.m function for solving a convection diffusion equation

% - grad div u + [1 0] grad u
% exact solution given (see documenatation)
% dirichlet boundary
% computational domain: square [0 1]^2

% define constants

NREFS = 4;            % number of refinement steps
a = 100;          % diffusion constant
plot = 0;             % plotting flag
output = 0;           % output flag

% solve problem

[errs, mw, msh, u] =exp_3(a,NREFS,plot,output);

% plot the solution

plot_Mesh(msh,'as')
plot_LFEfancy(u,msh);
plot_LFE(u,msh); colorbar;

disp(['MeshWidth = ',num2str(mw)]);

% clear memory

clear all;

