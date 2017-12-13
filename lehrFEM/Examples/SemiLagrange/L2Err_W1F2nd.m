function err = L2Err_W1F2nd(Mesh,u,QuadRule,FHandle,varargin)
% L2ERR_W1F Discretization error in L2 norm for W1F finite element.
%
%   ERR = L2ERR_W1F(MESH,U,QUADRULE,FHANDLE) computes the discretization
%   error between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
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
%   ERR = L2ERR_W1F(MESH,U,QUADRULE,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2Err_W1F(Mesh,u,QuadRule,F_Handle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  
  % Precompute shape function values at the quedrature points
  
  N = shap_W1F2nd(QuadRule.x);
    
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
    TK = transpose(inv(BK));

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
      
    u_EX = FHandle(x,varargin{:});
    u_FE = ( u(eidx(1)) * N(:,1:2) * p1 + ...
             u(eidx(2)) * N(:,3:4) * p2 + ...
             u(eidx(3)) * N(:,5:6) * p3+ ...
             u(eidx(1)+nEdges) * N(:,7:8) + ...
             u(eidx(2)+nEdges) * N(:,9:10) + ...
             u(eidx(3)+nEdges) * N(:,11:12))*TK;
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*sum(abs(u_EX-u_FE).^2,2))*det_BK;
      
  end
  
  err = sqrt(err);
  
return