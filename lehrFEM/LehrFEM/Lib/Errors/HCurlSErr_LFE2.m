function err = HCurlSErr_LFE2(Mesh,u,QuadRule,CURL_Handle,varargin)
% HCURLSERR_LFE2 Discretization error in S1 norm for LFE2 finite element.
%
%   ERR = HCURLSERR_LFE2(MESH,U,QUADRULE,CURL_HANDLE) computes the 
%   discretization error between the exact solution given by the function 
%   handle CURL_HANDLE and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = HCURLSERR_LFE2(MESH,U,QUADRULE,CURL_HANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution CURL_HANDLE.
%
%   Example:
%
%   err = HCurlSErr_LFE2(Mesh,u,QuadRule,CURL_Handle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  
  % Compute discretization error

  err = 0;
  for i = 1:nElements
       
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    vidxNeo = vidx + nCoordinates;
          
    % Compute element mapping

    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    bK = P1;
    BK = [P2-bK;P3-bK];
    det_BK = abs(det(BK));

    % Transform quadrature points

    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX = CURL_Handle(x,varargin{:} );
    u_FE = (u(vidx(1))*(P3(1)-P2(1)) + u(vidxNeo(1))*(P3(2)-P2(2)) + ...
            u(vidx(2))*(P1(1)-P3(1)) + u(vidxNeo(2))*(P1(2)-P3(2)) + ...
            u(vidx(3))*(P2(1)-P1(1)) + u(vidxNeo(3))*(P2(2)-P1(2)))/det_BK;
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*abs(u_EX-u_FE).^2)*det_BK;
      
  end
  
  err = sqrt(err);
  
return