function err = L2Norm_W1F2ndOutCircle(Mesh,u,QuadRule,R,C,varargin)
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
    r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
    xout =(r >= R);   
    
    u_FE = [xout, xout].*( u(eidx(1)) * N(:,1:2) * p1 + ...
             u(eidx(2)) * N(:,3:4) * p2 + ...
             u(eidx(3)) * N(:,5:6) * p3+ ...
             u(eidx(1)+nEdges) * N(:,7:8) + ...
             u(eidx(2)+nEdges) * N(:,9:10) + ...
             u(eidx(3)+nEdges) * N(:,11:12))*TK;
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*sum((u_FE).^2,2))*det_BK;
      
  end
  
  err = sqrt(err);
  
return