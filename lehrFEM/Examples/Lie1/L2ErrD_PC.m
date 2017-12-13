function err = L2ErrD_PC(Mesh,u,QuadRule,FHandle,varargin)
% L2ERRD_PC Discretization error in L2 norm for pecewise constant finite elements.
%
%   ERR = L2ERRD_LFE(MESH,U,QUADRULE,FHANDLE) computes the discretization
%   error between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = L2ERR_LFE(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2Err_LFE(Mesh,u,QuadRule,fhandle);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % barycentric refinement 

  b_Mesh=refine_BAR(Mesh);
  b_Mesh = add_Edge2Elem(b_Mesh);
  b_Mesh = add_Patches(b_Mesh);
  b_nElements=size(b_Mesh,1);
  b_U=zeros(b_nElements,1);

% representation of U on finer mesh  
  
  nCoordinates=size(Mesh.Coordinates,1);
  
  for i=1:nCoordinates
      b_AdjElements=b_Mesh.AdjElements(i,:);
      b_AdjElements=setdiff(b_AdjElements,0);
      b_nAdjElements=size(b_AdjElements,2);
      b_U(b_AdjElements)=u(i);
  end
 % error callculated on finer mesh
 
  err=L2Err_PC(b_Mesh,b_U,QuadRule,FHandle);
  
return