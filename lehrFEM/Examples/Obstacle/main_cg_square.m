% MAIN_CG_SQUARE
% Error analysis for the conjugate gradient method on a square mesh when
% solving the obstackle problem

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc
% Initializing constants

Shape = 'Sqr';      % Type of mesh
TOL = 1e-10;        % Stopping criterion 
MAXIT = inf;      % Maximum number of iterations
NREFS = 5;          % Number of red refinements

main_cg(Shape,NREFS,TOL,MAXIT);