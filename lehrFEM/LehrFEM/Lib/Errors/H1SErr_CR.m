function err = H1SErr_CR(Mesh,u,QuadRule,FHandle,varargin)
% H1SERR_CR Discretization error in H1 semi-norm for quadratic finite element.
%
%   ERR = H1SERR_CR(MESH,U,QUADRULE,FHANDLE) computes the discretization
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
%   ERR = H1SERR_CR(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = H1SErr_CR(Mesh,u,QuadRule,F_Handle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
    
  % Precompute gradients of shape functions
  
  grad_N = grad_shap_CR(QuadRule.x);
  
  % Compute discretization error 
   
  err = 0;
  eidx = zeros(1,3);
  for i= 1:nElements
        
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(2),vidx(3));
    eidx(2) = Mesh.Vert2Edge(vidx(3),vidx(1));
    eidx(3) = Mesh.Vert2Edge(vidx(1),vidx(2));
    
    % Compute element mapping
      
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    
    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
        
    % Evaluate solutions
      
    grad_u_EX = FHandle(x,varargin{:});
    grad_u_FE = (u(eidx(1))*grad_N(:,1:2) + ...
                 u(eidx(2))*grad_N(:,3:4) + ...
                 u(eidx(3))*grad_N(:,5:6))*transpose(inv_BK);
      
    % Compute error on the current element
     
    err = err+sum(QuadRule.w.*sum(abs(grad_u_FE-grad_u_EX).^2,2))*det_BK;
      
  end
  
  err = sqrt(err);
    
return

