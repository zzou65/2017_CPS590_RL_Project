function err = HveloCurlSErr_W1F(Mesh,u,QuadRule,VCURL_Handle,V_Handle,varargin)
% HCURLSERR_W1F Discretization error in S1 norm for W1F finite element.
%
%   ERR = HCURLSERR_W1F(MESH,U,QUADRULE,CURL_HANDLE) computes the 
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
%   ERR = HCURLSERR_W1F(MESH,U,QUADRULE,CURL_HANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution CURL_HANDLE.
%
%   Example:
%
%   err = HCurlSErr_W1F(Mesh,u,QuadRule,CURL_Handle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Compute discretization error

  err = 0;
  eidx = zeros(1,3);
  for i = 1:nElements
       
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(2),vidx(3));
    eidx(2) = Mesh.Vert2Edge(vidx(3),vidx(1));
    eidx(3) = Mesh.Vert2Edge(vidx(1),vidx(2));
          
    % Compute element mapping

    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
        Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));

    % Determine the orientation

    if(Mesh.Edges(eidx(1),1) == vidx(2))
        p1 = 1;
    else
        p1 = -1;
    end

    if(Mesh.Edges(eidx(2),1) == vidx(3))
        p2 = 1;
    else
        p2 = -1;
    end

    if(Mesh.Edges(eidx(3),1) == vidx(1))
        p3 = 1;
    else
        p3 = -1;
    end
        
    % Transform quadrature points

    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
    Fval = V_Handle(x);
    Fval = [Fval(:,2) -Fval(:,1)];
    u_EX = VCURL_Handle(x,varargin{:});
    u_FE = -2/det_BK*(u(eidx(1))*p1+u(eidx(2))*p2+u(eidx(3))*p3)*Fval;
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*sum((u_EX-u_FE).^2,2))*det_BK;
      
  end
  
  err = sqrt(err);
  
return