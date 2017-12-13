% Test routine for 1D hpDG discretization schemes.

% Copyright 2007-2007 Patrick Meury & Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize mesh
  
  X0 = -1;                  % Left end point of interval
  X1 =  1;                  % Right end point of interval
  NPOINTS = 20;             % Number of grid points
  PMAX = 0;                 % Maximal polynomial degree
  F = @(x,varargin)exp(x);  % Interpolated function
  
  % Initialize mesh
  
  Coordinates = X0 + (X1-X0)/(NPOINTS-1)*(0:(NPOINTS-1));
  
  % Initialize polynomial degrees
  
  p = PMAX*ones(1,NPOINTS-1);
  
  % Initialize quadrature rule and shape functions
  
  qr = gauleg(-1,1,200);
  Shap = shap_Leg_1D(qr.x,PMAX);
  
  % Assemble mass matrices and load vector
  
  M = assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,qr,Shap);
  L = assemLoad_hpDG_1D(Coordinates,p,qr,Shap,F);

  % Solve the linear system
  
  u = M\L;
  
  % Plot solution
  
  plot_hpDG_1D(Coordinates,p,u,@shap_Leg_1D);
  
  % Clear memory
  
  clear all;