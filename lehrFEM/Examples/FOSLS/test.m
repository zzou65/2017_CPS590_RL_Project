% Run script for piecewise linear finite element solver.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
%   NREFS = 3;               % Number of red refinement steps
%   F_HANDLE = @f_LShap;     % Right hand side source term
%   GD_HANDLE = @g_D_LShap;  % Dirichlet boundary data
%   %GN_HANDLE = @g_N_LShap;  % Neumann boundary data
%  
%   % Initialize mesh
%   
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
%   
 clear Mesh;
    
    % Initialize constant
    
    NREFS = 4;
    F_Handle = @(x,varargin)(pi^2+1)*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
%     F_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];
%     GD_Handle = @(x,varargin)[-x(:,2) x(:,1)];

    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
    
        Mesh = refine_REG(Mesh);
    
    [IM,JM,M] = assemMat_W1F(Mesh,@MASS_W1F,0,P7O6());
    A = sparse([IM],[JM],[M]);
    L = assemLoad_W1F(Mesh,P7O6(),F_Handle);
    
    U= A\L;
    L2Err_TR(NewMesh,U,P7O6(),F_Handle);
    end
    
  
  clear all;
  