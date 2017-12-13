% Computes solutions to the steady Stokes problem on the unit square using
% piecewise quadratic finite elements for the velocity and piecewise constants
% for the pressure.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 5;             % Number of red refinements
  
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

  plot_Mesh(Mesh,'as')

  Loc = get_BdEdges(Mesh);
  Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
  FixedPos = zeros(size(Mesh.Coordinates,1),1);
  FixedPos(Loc) = 1;
  Mesh = jiggle(Mesh,FixedPos);

  plot_Mesh(Mesh,'as')
  
  for NU = .01:.01:.09

      GD = @(x,varargin)[cos(pi/2*(x(:,1)+x(:,2))) -cos(pi/2*(x(:,1)+x(:,2)))];                  % Dirichlet boundary data
      F = @(x,varargin)[ NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))-pi/2*cos(pi/2*(x(:,1)-x(:,2))) ...  % Right hand side source
          -NU*pi^2/2*cos(pi/2*(x(:,1)+x(:,2)))+pi/2*cos(pi/2*(x(:,1)-x(:,2)))];

      % Assemble stiffness matrix and load vector

      A = assemMat_Stokes_P1P0(Mesh,@STIMA_Stokes_P1P0,NU,P7O6());
      L = assemLoad_Stokes_P1P0(Mesh,P7O6(),F);

      % Incorporate Dirichlet boundary data

      [U,FreeDofs] = assemDir_Stokes_P1P0(Mesh,-1,GD);
      L = L - A*U;

      % Solve the linear system

      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);

  end
   