function err = H1Err_PBD(Mesh,u,QuadRule,FHandle,varargin)
% H1ERR_PBD Discretization error with respect to H1 norm for quadratic
%           finite elements with parabolic boundary approximation.
%
%   ERR = H1ERR_PBD(MESH,U,QUADRULE,FHANDLE) computes the discretization
%   error between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of
%                 elements.
%    DELTA        P-by-1 matrix specifying the boundary correction term on
%                 every edge.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = H1ERR_PBD(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = H1SErr_PBD(Mesh,u,QuadRule,FHandle);
%
%   See also shap_QFE, grad_shap_QFE.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  Rot = [0 -1; 1 0];
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Precompute shape functions
 
  N = shap_QFE(QuadRule.x);
  grad_N = grad_shap_QFE(QuadRule.x);
      
  % Compute discretization error
  
  err = 0;
  eidx = zeros(1,3);
  for i = 1:nElements
       
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(1),vidx(2)) + nCoordinates;
    eidx(2) = Mesh.Vert2Edge(vidx(2),vidx(3)) + nCoordinates;
    eidx(3) = Mesh.Vert2Edge(vidx(3),vidx(1)) + nCoordinates;
    
    % Compute standard element mapping  
        
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    
    delta = Mesh.Delta(Mesh.Vert2Edge(vidx(2),vidx(3)));
    if(delta > eps)
             
      % Compute curved element mapping
      
      normal = Mesh.Coordinates(vidx(3),:)-Mesh.Coordinates(vidx(2),:);
      normal = normal*Rot/norm(normal);
      
      x = QuadRule.x*BK + ones(nPts,1)*bK + 4*delta*N(:,1).*N(:,2)*normal;
      
      % Evaluate exact solution
      
      [u_EX,grad_u_EX] = FHandle(x,varargin{:});    
      
      % Compute L2 discretization error
      
      for j = 1:nPts
        Jac = transpose(BK) + 4*delta*transpose(normal) * ...
              [-QuadRule.x(j,2) 1-QuadRule.x(j,1)-2*QuadRule.x(j,2)];
        inv_Jac = inv(Jac);
        det_Jac = abs(det(Jac));
        u_FE = u(vidx(1))*N(j,1)+u(vidx(2))*N(j,2)+u(vidx(3))*N(j,3) + ...
               u(eidx(1))*N(j,4)+u(eidx(2))*N(j,5)+u(eidx(3))*N(j,6);
        grad_u_FE = (u(vidx(1))*grad_N(j,1:2) + ...
                     u(vidx(2))*grad_N(j,3:4) + ...
                     u(vidx(3))*grad_N(j,5:6) + ... 
                     u(eidx(1))*grad_N(j,7:8) + ...
                     u(eidx(2))*grad_N(j,9:10) + ...
                     u(eidx(3))*grad_N(j,11:12))*inv_Jac;
        err = err + QuadRule.w(j)*(abs(u_EX(j)-u_FE)^2 + ...
                                   norm(grad_u_EX(j,:)-grad_u_FE))^2*det_Jac;
      end    
      
    else
    
      % Compute straight element mapping
    
      inv_BK = inv(BK);
      det_BK = abs(det(BK));
      x = QuadRule.x*BK+ones(nPts,1)*bK;
        
      % Evaluate solutions
      
      [u_EX,grad_u_EX] = FHandle(x,varargin{:});
      u_FE = u(vidx(1))*N(:,1)+u(vidx(2))*N(:,2)+u(vidx(3))*N(:,3) + ...
             u(eidx(1))*N(:,4)+u(eidx(2))*N(:,5)+u(eidx(3))*N(:,6);
      grad_u_FE = (u(vidx(1))*grad_N(:,1:2) + ...
                   u(vidx(2))*grad_N(:,3:4) + ...
                   u(vidx(3))*grad_N(:,5:6) + ...
                   u(eidx(1))*grad_N(:,7:8) + ...
                   u(eidx(2))*grad_N(:,9:10) + ...
                   u(eidx(3))*grad_N(:,11:12))*transpose(inv_BK);
    
      % Compute L2 discretization error
      
      err = err+sum(QuadRule.w.*(abs(u_EX-u_FE).^2 + ...
                                 sum(abs(grad_u_EX-grad_u_FE).^2,2)))*det_BK;
      
    end
  end
  err = sqrt(err);
  
return