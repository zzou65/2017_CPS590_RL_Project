% Computes solutions to the steady Stokes problem on the unit square using
% Crouzeix-Raviart finite elements for the velocity and piecewise constant
% finite elements for the pressure.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 2;                                                                                     % Number of red refinements
  NU = 1;                                                                                        % Viscosity
  GD_HANDLE = @(x,varargin)[cos(pi/2*(x(:,1)+x(:,2))) -cos(pi/2*(x(:,1)+x(:,2)))];               % Dirichlet boundary data
  F_HANDLE = @(x,varargin)[ NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))-pi/2*cos(pi/2*(x(:,1)-x(:,2))) ...  % Right hand side source
                           -NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))+pi/2*cos(pi/2*(x(:,1)-x(:,2)))];
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS
   Mesh = refine_REG(Mesh);     
  end
  
  % Assemble stiffness matrix and load vector
 
  A = assemMat_Stokes_CRP0(Mesh,@STIMA_Stokes_CRP0,NU,P7O6());
  L = assemLoad_Stokes_CRP0(Mesh,P7O6(),F_HANDLE);
  
  % Incorporate Dirichlet boundary data
 
  [U,FreeDofs] = assemDir_Stokes_CRP0(Mesh,-1,GD_HANDLE);
  L = L - A*U; 
  
  % Solve the linear system
 
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
  % Plot out solution
    
  plot_Stokes(U,Mesh,'CRP0');
  title('{\bf Steady Stokes equation (Crouzeix-Raviart elements)}');
  xlabel(['{\bf # Dofs  :  ' int2str(size(U,1)) '}']);
  colorbar;
  
  % Output .eps files
  
  print('-depsc','func_CRP0.eps')
  
  % Clear memory
  
%   clear all;
  