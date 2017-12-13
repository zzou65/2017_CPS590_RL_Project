function J = Flux_QFE(u,BdFlags,Mesh,QuadRule,SigmaHandle,varargin)
% FLUX_LFE Flux throu specified part of the boundary.
%
%   J = FLUX_LFE(U,BDFLAGS,MESH,QUADRULE,SIGMAHANDLE) computes the value
%   of the flux throu the boundary specified by the integers BDFLAGS of
%   the finite element solution U on the struct MESH using the quadrature
%   rule QUADRULE.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.   
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    EDGE2ELEM   P-by-2 matrix connecting edges to elements. The first
%                column specifies the left hand side element where the
%                second column specifies the right hand side element. 
%    EDGELOC     P-by-3 matrix connecting egdes to local edges of
%                elements.
%
%   J = FLUX_LFE(U,BDFLAGS,MESH,QUADRULE,SIGMAHANDLE,SIGMAPARAM) also
%   handles the variable length argument list SIGMAPARAM to the function
%   handle SIGMAPARAM during assembly.
%
%   Example:
%
%   SigmaHandle = @(x,varargin)ones(size(x,1),1);
%   J = Flux_LFE(u,-1,Mesh,QuadRule,SigmaHandle);
%
%   See also GET_BDEDGES, SHAP_LFE.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nGauss = size(QuadRule.w,1);
  Rot = [0 -1; 1 0];
  
  % Precompute shape functions
  
  grad_N = grad_shap_QFE([QuadRule.x zeros(nGauss,1)]);
  
  J = 0;
  for j1 = BdFlags
  
    % Extract Neumann edges
  
    Loc = get_BdEdges(Mesh);           
    Loc = Loc(Mesh.BdFlags(Loc) == j1);
     
    for j2 = Loc'
        
      % Compute element map
      
      if(Mesh.Edge2Elem(j2,1))
        
        % Match orientation to left hand side element  
        
        Elem = Mesh.Edge2Elem(j2,1);
        id_s = Mesh.Elements(Elem,rem(Mesh.EdgeLoc(j2,1),3)+1);
        id_e = Mesh.Elements(Elem,rem(Mesh.EdgeLoc(j2,1)+1,3)+1);
        
      else
          
        % Match orientation to right hand side element  
          
        Elem = Mesh.Edge2Elem(j2,2);
        id_s = Mesh.Elements(Elem,rem(Mesh.EdgeLoc(j2,2),3)+1);
        id_e = Mesh.Elements(Elem,rem(Mesh.EdgeLoc(j2,2)+1,3)+1);
        
      end
      Q0 = Mesh.Coordinates(id_s,:);
      Q1 = Mesh.Coordinates(id_e,:);
      x = ones(nGauss,1)*Q0+QuadRule.x*(Q1-Q0);
      dS = norm(Q1-Q0);
      
      % Compute outer unit normals
      
      normal = ones(nGauss,1)*((Q1-Q0)/dS*Rot);
     
      % Compute function values on the current edge
      
      SigmaVal = SigmaHandle(x,Mesh.ElemFlag(Elem),varargin{:});
      vidx = Mesh.Elements(Elem,:);    
      idx = [vidx ...
             Mesh.Vert2Edge(vidx(1),vidx(2))+nCoordinates ...
             Mesh.Vert2Edge(vidx(2),vidx(3))+nCoordinates ...
             Mesh.Vert2Edge(vidx(3),vidx(1))+nCoordinates];
      BK = [Mesh.Coordinates(vidx(2),:)-Mesh.Coordinates(vidx(1),:); ...
            Mesh.Coordinates(vidx(3),:)-Mesh.Coordinates(vidx(1),:)];
      inv_BK_t = transpose(inv(BK));
      grad_U = (u(idx(1))*grad_N(:,1:2)+u(idx(2))*grad_N(:,3:4) + ...
                u(idx(3))*grad_N(:,5:6)+u(idx(4))*grad_N(:,7:8) + ...
                u(idx(5))*grad_N(:,9:10)+u(idx(6))*grad_N(:,11:12))*inv_BK_t;
      
      % Update flux      
            
      J = J + sum(QuadRule.w.*SigmaVal.*sum(grad_U.*normal,2))*dS;
      
    end    
  end

return