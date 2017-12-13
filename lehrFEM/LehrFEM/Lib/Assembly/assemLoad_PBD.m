function L = assemLoad_PBD(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_PBD Assemble quadratic FE contributions.
%
%   L = ASSEMLOAD_PBD(MESH,QUADRULE,FHANDLE) assembles the global load
%   vector for the load data given by the function handle EHANDLE and
%   parabolic boundary corrections.
%
%   The struct MESH must at least contain the following fields:
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
%   L = ASSEMLOAD_PBD(COORDINATES,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemLoad_PBD(Mesh,P7O6(),FHandle);
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  Rot = [0 -1; 1 0];
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates+nEdges,1);
  
  % Precompute shape functions
  
  N = shap_QFE(QuadRule.x);
  
  eidx = zeros(1,3);
  lam_1 = 1-QuadRule.x(:,1)-QuadRule.x(:,2);
  lam_2 = QuadRule.x(:,2);
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
      x = QuadRule.x*BK + ones(nPts,1)*bK + 4*delta*lam_1.*lam_2*normal;
      
      % Compute load data
      
      FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
      
      % Add curved contributins to global load vector
      
      for i = 1:nPts
        Jac = transpose(BK) + 4*delta*transpose(normal) * ...
              [-QuadRule.x(i,2) 1-QuadRule.x(i,1)-2*QuadRule.x(i,2)];  
        det_Jac = abs(det(Jac));
        L(vidx(1)) = L(vidx(1)) + QuadRule.w(i)*FVal(i)*N(i,1)*det_Jac;
        L(vidx(2)) = L(vidx(2)) + QuadRule.w(i)*FVal(i)*N(i,2)*det_Jac;
        L(vidx(3)) = L(vidx(3)) + QuadRule.w(i)*FVal(i)*N(i,3)*det_Jac;
        L(eidx(1)) = L(eidx(1)) + QuadRule.w(i)*FVal(i)*N(i,4)*det_Jac;
        L(eidx(2)) = L(eidx(2)) + QuadRule.w(i)*FVal(i)*N(i,5)*det_Jac;
        L(eidx(3)) = L(eidx(3)) + QuadRule.w(i)*FVal(i)*N(i,6)*det_Jac;
      end    

    else
    
      % Compute straight element mapping  
      
      det_BK = abs(det(BK));
      x = QuadRule.x*BK + ones(nPts,1)*bK;
      
      % Compute load data
      
      FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
      
      % Add straight contributins to global load vector
      
      L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1))*det_BK;
      L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2))*det_BK;
      L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3))*det_BK;     
      L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*FVal.*N(:,4))*det_BK;
      L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*FVal.*N(:,5))*det_BK;
      L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*FVal.*N(:,6))*det_BK;
      
    end  
  end
  
return