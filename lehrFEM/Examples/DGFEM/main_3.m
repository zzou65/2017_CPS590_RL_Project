% Run script for discontinuous Galerkin finite element solver

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   clear all
  NREFS = 6;                                % Number of red refinement steps
  G = @(x,varargin)-4*ones(size(x,1),1);    % Right hand side load data
  UD = @(x,varargin)x(:,1).^2+x(:,2).^2;    % Dirichlet boundary data
  S = 1;                                    % Symmetric (+1) or antisymmetric (-1) discretization
  SIGMA = @(P0,P1,varargin)10/norm(P1-P0);  % Edge weight function
  
  % Initialize mesh
  
  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
    
  % Assemble matrices and load vectors (discontinuous Lagrangian elements)
   
  QuadRule_1D = gauleg(0,1,2);
  QuadRule_2D = P3O3();
  
  [I1,J1,Avol] = assemMat_Vol_DG(Mesh,@STIMA_Lapl_Vol_DGLFE);
  
  [I2,J2,Jinn] = assemMat_Inn_DG(Mesh,@STIMA_InnPen_DGLFE,SIGMA);
  [I2,J2,Ainn] = assemMat_Inn_DG(Mesh,@STIMA_Inn_DGLFE,S);
    
  [I3,J3,Jbnd] = assemMat_Bnd_DG(Mesh,@STIMA_BndPen_DGLFE,SIGMA);  
  [I3,J3,Abnd] = assemMat_Bnd_DG(Mesh,@STIMA_Bnd_DGLFE,S);
  
  Lvol = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGLFE,QuadRule_2D,G);
  Lbnd = assemLoad_Bnd_DG(Mesh,@LOAD_Bnd_DGLFE,QuadRule_1D,S,SIGMA,UD);
    
  % Create system matrix
       
  A = sparse([J1; J2; J3; J2; J3], ...
             [I1; I2; I3; I2; I3], ...
             [Avol; Ainn; Abnd; Jinn; Jbnd]);
  L = Lvol + Lbnd;
         
  % Solve the linear system
  
  U = A\L;
  plot_DGLFE(U,Mesh);
  colorbar;
      
  % Clear memory
  
  clear all;
  