 % Run script for piecewise linear finite element solver.

 %   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
 %   SAM - Seminar for Applied Mathematics
 %   ETH-Zentrum
 %   CH-8092 Zurich, Switzerland
  
  NREFS =4; 
  F_HANDLE = @f_;        % Right hand side source term
 % gn_HANDLE = @gn_;     % neumann data
  g_HANDLE = @g_;        % impedance data 
  DHANDLE = @dist_circ;  % Signed distance function
  C = [0 0];             % Center of the circle
  R = 1;                 % Radius of the circle
  sol= @(x,varargin)(1+x(:,1).^2+x(:,2).^2).^(-1/2);  

  % Load mesh from file

  Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
  Mesh = shift(Mesh,C);

  % Add edge data structure

  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;

  % Do NREFS uniform refinement steps
  plot_Mesh(Mesh,'petas');
  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
  %  plot_Mesh(Mesh,'as');
  end  
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1); 
  % Assemble stiffness matrix and load vector
 
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  prec=1e-6;
  points=10;
  B = assemMat_Bnd_Robin(Mesh,g_HANDLE, prec, points);
  A=A-B;
   % Incorporate Impedance type adds to system
  L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE,1);
   
 

  % Incorporate Neumann boundary data
  
  % L = assemNeu_LFE(Mesh,-1,L,gauleg(0,1,4),GN_HANDLE);
  
 
  % Solve the linear system
 
  U = A\L;
    
  % Plot out solution
    
  plot_LFE(U,Mesh);
  colorbar;
  plot_LFE(sol(Mesh.Coordinates()),Mesh);
  colorbar;
  
  % Clear memory
  
  clear all;
  