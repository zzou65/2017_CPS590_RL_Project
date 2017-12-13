% MAIN_PSOR_SQUARE
% Error analysis for the projected sor method on a square mesh when
% solving the obstackle problem

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc
% Initializing constants

Shape = 'Sqr';          % Type of mesh
TOL = 1e-10;            % Stopping criterion 
MAXIT = 10000;          % Maximum number of iterations
NREFS = 5;              % Number of red refinements
omega = 2/(1+sqrt(2));  % SOR extrapolation factor

main_psor(Shape,NREFS,TOL,MAXIT,omega);