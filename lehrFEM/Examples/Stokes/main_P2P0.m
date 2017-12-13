% LehrFEm driver script for computing solutions of the steady Stokes problem on the unit square
% using piecewise quadratic finite elements for the velocity and piecewise constants for the
% pressure.
  NREFS = 4; % Number of red refinements
  NU = 1;    % Viscosity
  % Dirichlet boundary data
  GD_HANDLE = @(x,varargin)[cos(pi/2*(x(:,1)+x(:,2))) -cos(pi/2*(x(:,1)+x(:,2)))];
  % Right hand side source (force field)
  F_HANDLE = @(x,varargin)[ sin(pi*x(:,1)) zeros(size(x(:,1))) ];
  
  % Initialize mesh
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS, Mesh = refine_REG(Mesh); end
  
  % Assemble Galerkin matrix and load vector
  A = assemMat_Stokes_P2P0(Mesh,@STIMA_Stokes_P2P0,NU,P7O6());
  L = assemLoad_Stokes_P2P0(Mesh,P7O6(),F_HANDLE);
  
  % Incorporate Dirichlet boundary data
  [U,FreeDofs] = assemDir_Stokes_P2P0(Mesh,-1,GD_HANDLE); L = L - A*U;
  
  % Solve the linear system (direct solver)
  U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
  % Plot and print solution
  plot_Stokes(U,Mesh,'P2P0');
  title('{\bf Steady Stokes equation (P2 elements)}');
  xlabel(['{\bf # Dofs  :  ' int2str(size(U,1)) '}']);
  colorbar;
  print('-depsc','func_P2P0.eps')
