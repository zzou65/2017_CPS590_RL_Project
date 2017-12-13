function err = Helicity2D_W1F(Mesh,u,QuadRule,V_Handle,varargin)

%   Copyright 2010 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  
  % Compute discretization error
  
  % Precompute shape function values at the quedrature points
  
  N = shap_W1F(QuadRule.x);

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
    velo = V_Handle(x,varargin{:});
    velo = [velo(:,2), -velo(:,1)];
      
    u_FE = ( u(eidx(1)) * N(:,1:2) * p1 + ...
        u(eidx(2)) * N(:,3:4) * p2 + ...
        u(eidx(3)) * N(:,5:6) * p3)*TK;
    
    
    du_FE = -2/det_BK*(u(eidx(1))*p1+u(eidx(2))*p2+u(eidx(3))*p3)*ones(nPts,1);
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*du_FE.*sum(velo.*u_FE,2))*det_BK;
      
  end
  
  err = err;
  
return