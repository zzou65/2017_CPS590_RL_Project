% Run script for 1D finite element solver.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum

  % Initialize constants

  XL = 0;                                       % Left hand-side boundary
  XR = 1;                                       % Right hand-side boundary
  MESH_WIDTH = 2.^(-2:-1:-10);                  % Mesh resolution 
  F = @(x,varargin)8*ones(size(x));             % Right hand-side source term
  GD = @(x,varargin)zeros(size(x));             % Dirichlet boundary data
  GN = @(x,varargin)4*ones(size(x));            % Neumann boundary data
  UEX_1 = @(x,varargin)4*x.*(1-x);              % Exact solution for L2 norm   
  UEX_2 = @(x,varargin)deal(4*x.*(1-x),4-8*x);  % Exact solution for H1 norm
  
  % Compute discretization errors
  
  h = zeros(size(MESH_WIDTH));
  L2_err = zeros(size(MESH_WIDTH));
  H1_err = zeros(size(MESH_WIDTH));
  for i = 1:size(MESH_WIDTH,2)
        
    % Initialize mesh
  
    Coordinates = transpose(XL:MESH_WIDTH(i):XR);
  
    % Compute stiffness matrix and load vector
  
    A = assemMat_P1_1D(Coordinates,@STIMA_Lapl_P1_1D);
    L = assemLoad_P1_1D(Coordinates,gauleg(0,1,2),F);
  
    % Incorporate Neumann boundary data
 
    L = assemNeu_P1_1D(Coordinates,1,L,GN);
  
    % Incorporate Dirichlet boundary data
 
    [U,FreeDofs] = assemDir_P1_1D(Coordinates,size(Coordinates,1),GD);
    L = L - A*U;
    
    % Solve the linear system
  
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % Compute discretiyation error
    
    L2_err(i) = L2Err_P1_1D(Coordinates,U,gauleg(0,1,2),UEX_1);
    H1_err(i) = H1Err_P1_1D(Coordinates,U,gauleg(0,1,2),UEX_2);
    
  end
    
  % Generate figures
  
  fig = figure('Name','Discretization errors');
  plot(MESH_WIDTH,L2_err,'r-', ...
       MESH_WIDTH,H1_err,'b-', ...
       MESH_WIDTH,L2_err,'k+', ...
       MESH_WIDTH,H1_err,'k+');
  grid('on');
  set(gca,'XScale','log','XDir','reverse','YScale','log');
  title('{\bf Discretization error for linear finite elements}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('L^2 norm','H^1 norm','Location','SouthEast');
  p = polyfit(log(MESH_WIDTH),log(L2_err),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(MESH_WIDTH),log(H1_err),1);
  add_Slope(gca,'NorthWest',p(1));
   
  % Clear memory
  
  clear all;
    