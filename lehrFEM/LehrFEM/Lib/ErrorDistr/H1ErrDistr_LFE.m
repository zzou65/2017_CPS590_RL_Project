function err = H1ErrDistr_LFE(Mesh,u,QuadRule,FHandle,varargin)
% H1ERRDISTR_LFE Discretization error distribution in H1 norm for linear 
% finite element
%
%   ERR = H1ERRDISTR_LFE(MESH,U,QUADRULE,FHANDLE) computes the discretization
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
%   ERR = H1ERRDISTR_LFE(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = H1ErrDistr_LFE(Mesh,u,QuadRule,fhandle);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  % Inititialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  err = zeros(nElements,1);
  
  % Precompute gradients and values of shape functions
  
  N = shap_LFE(QuadRule.x);
  grad_N = grad_shap_LFE(QuadRule.x);
    
  % Compute discretization error
 
  for i= 1:nElements
          
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
      
    % Compute element mapping
      
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
      
    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
        
    % Evaluate solutions
      
    [u_EX,grad_u_EX] = FHandle(x,varargin{:});
    u_FE = u(vidx(1))*N(:,1)+u(vidx(2))*N(:,2)+u(vidx(3))*N(:,3);
    grad_u_FE = (u(vidx(1))*grad_N(:,1:2)+ ...
                 u(vidx(2))*grad_N(:,3:4)+ ...
                 u(vidx(3))*grad_N(:,5:6))*transpose(inv_BK);
      
    % Compute error on the current element
      
    err(i) = sum(QuadRule.w.*(abs(u_EX-u_FE).^2+...
                    sum(abs(grad_u_FE-grad_u_EX).^2,2)))*det_BK;
      
  end 
  
  err = sqrt(err);
    
return