% Discretization error of the mean value for piecewise quadratic and linear
% finite elements for the Poisson equation with zero Dirichlet boundary
% conditions on the unit square.
 
%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 7;                                                      % Number of red refinement steps
  F_HANDLE = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));  % Right hand side source term
  GD_HANDLE = @(x,varargin)zeros(size(x,1),1);                    % Dirichlet boundary data
  FEX = 4/pi^2;                                                   % Exact value of functional
  JIG = 1;                                                        % Jiggle parameter
    
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Compute discretization error on a series of meshes
    
  h = zeros(1,NREFS);
  N_LFE = zeros(1,NREFS);
  N_QFE = zeros(1,NREFS);
  MeanErr_LFE = zeros(1,NREFS);
  MeanErr_QFE = zeros(1,NREFS);
  for i = 1:NREFS
      
    % Do red mesh refinement  
      
    Mesh = refine_REG(Mesh);    
    Mesh = add_Edge2Elem(Mesh);
    
    % Mesh preprocessing
    
    switch (JIG)
      case 1
        NewMesh = Mesh;
      case 2
        Loc = get_BdEdges(Mesh);
        Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        NewMesh = jiggle(Mesh,FixedPos);
      case 3
        Loc = get_BdEdges(Mesh);
        Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        NewMesh = smooth(Mesh,FixedPos);
    end
    
    % Assemble stiffness matrix, load vector and incorporate BC
  
    A_QFE = assemMat_QFE(NewMesh,@STIMA_Lapl_QFE);
    L_QFE = assemLoad_QFE(NewMesh,P7O6(),F_HANDLE);
    
    A_LFE = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
    L_LFE = assemLoad_LFE(NewMesh,P7O6(),F_HANDLE);
         
    [U_LFE,FreeDofs_LFE] = assemDir_LFE(NewMesh,-1,GD_HANDLE);
    L_LFE = L_LFE - A_LFE*U_LFE;
    
    [U_QFE,FreeDofs_QFE] = assemDir_QFE(NewMesh,-1,GD_HANDLE);
    L_QFE = L_QFE - A_QFE*U_QFE;
    
    % Solve the linear system
  
    U_LFE(FreeDofs_LFE) = A_LFE(FreeDofs_LFE,FreeDofs_LFE)\L_LFE(FreeDofs_LFE);
    U_QFE(FreeDofs_QFE) = A_QFE(FreeDofs_QFE,FreeDofs_QFE)\L_QFE(FreeDofs_QFE);
            
    % Compute mean error
    
    MeanErr_LFE(i) = abs(FEX - mean_LFE(U_LFE,NewMesh));
    MeanErr_QFE(i) = abs(FEX - mean_QFE(U_QFE,NewMesh));
    
    % Compute mesh width and number dofs
    
    h(i) = get_MeshWidth(NewMesh);
    N_LFE(i) = size(NewMesh.Coordinates,1);
    N_QFE(i) = size(NewMesh.Coordinates,1)+size(NewMesh.Edges,1);
    
  end 
  
  % Plot out mean error against h mesh width and add slope triangles
  
  fig = figure('Name','Discretization errors');
  plot(h,MeanErr_QFE,'r-', ...
       h,MeanErr_LFE,'b-', ...
       h,MeanErr_QFE,'k+', ...
       h,MeanErr_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization error of mean value}');
  xlabel('{\bf Mesh width [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(h),log(MeanErr_LFE),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(h),log(MeanErr_QFE),1);
  add_Slope(gca,'South',p(1));
  
  % Plot out mean error against N number of dofs and add slope triangles
  
  fig = figure('Name','Discretization errors');
  plot(N_QFE,MeanErr_QFE,'r-', ...
       N_LFE,MeanErr_LFE,'b-', ...
       N_QFE,MeanErr_QFE,'k+', ...
       N_LFE,MeanErr_LFE,'k+');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Discretization error of mean value}');
  xlabel('{\bf Number of dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  
  legend('Quadratic FE','Linear FE','Location','NorthEast');
  p = polyfit(log(N_QFE),log(MeanErr_LFE),1);
  add_Slope(gca,'North',p(1));
  p = polyfit(log(N_LFE),log(MeanErr_QFE),1);
  add_Slope(gca,'South',p(1));
    
  % Clear memory
  
  clear all;
  