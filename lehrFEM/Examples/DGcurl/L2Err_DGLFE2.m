function err = L2Err_DGLFE2(Mesh,u,QuadRule,FHandle1,FHandle2,varargin)
% L2ERR_DGLFE2 Discretization error measured in the L2 norm for vectorial discontinuous
% galerkin finite elements.
%
%   ERR = L2ERR_DGLFE2(MESH,U,QUADRULE,FHANDLE1,FHANDLE2) computes the discretization
%   error between the exact solution given by the function handle FHANDLE1
%   FHANDLE2 and the finite element solution U on the struct MESH.
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
%   ERR = L2ERR_DGLFE(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
   % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  % Precompute shape function values at the quedrature points
  
  N = shap_DGLFE(QuadRule.x);
    
  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    idx = 6*(i-1)+[1 2 3 4 5 6];
          
    % Compute element mapping  
        
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));

    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX1 = FHandle1(x,ElemFlag(i),varargin{:});
    u_EX2 = FHandle2(x,ElemFlag(i),varargin{:});
    u_FE1 = u(idx(1))*N(:,1) + u(idx(3))*N(:,2) + u(idx(5))*N(:,3);
    u_FE2 = u(idx(2))*N(:,1) + u(idx(4))*N(:,2) + u(idx(6))*N(:,3);
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*abs(u_EX1-u_FE1).^2)*det_BK+sum(QuadRule.w.*abs(u_EX2-u_FE2).^2)*det_BK;
      
  end
  
  err = sqrt(err);
  
return