function err = L2Err_PC(Mesh,u,QuadRule,FHandle,varargin)
% L2ERR_PC Discretization error in L2 norm for pecewise constant finite elements.
%
%   ERR = L2ERR_PC(MESH,U,QUADRULE,FHANDLE) computes the discretization
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
%   ERR = L2ERR_PC(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2Err_LFE(Mesh,u,QuadRule,fhandle);

%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Precompute shape functions
  
  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract vertex numbers
    vidx=Mesh.Elements(i,:);
      
    % Compute element mapping  
        
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));

    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX = FHandle(x,varargin{:});
    u_FE = u(i)*ones(size(x,1),1);
      
    % Compute error on current element
    
    err = err+sum(QuadRule.w.*abs(u_EX-u_FE).^2)*det_BK;
      
  end
  err = sqrt(err);
  
return