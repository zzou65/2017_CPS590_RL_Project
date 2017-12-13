function L = assemNeu_hp(Mesh,Elem2Dof,L,BdFlags,QuadRule,Shap,FHandle,varargin)
% ASSEMNEU_HP Neumann boundary conditions.
%
%   L = ASSEMNEU_HP(MESH,ELEM2DOF,L,BDFLAGS,QUADRULE,SHAP,FHANDLE) adds the
%   Neumann boundary data with the data given by FHANDLE onto the load
%   vector L. 
%
%   The struct ELEM2DOF describes the element to dof mapping. The boundary
%   condition is only enforced on the edges whose boundary flag is equal to
%   the integer BDFLAG.
%
%   L = ASSEMNEU_HP(MESH,ELEM2DOF,L,BDFLAGS,QUADRULE,SHAP,FHANDLE,FPARAM)
%   also handles the variable length argument list FPARAM to the boundary
%   data function FHANDLE.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.   
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%    EDGE2ELEM   P-by-2 matrix connecting edges to elements. The first
%                column specifies the left hand side element where the
%                second column specifies the right hand side element.
%    EDGELOC     P-by-3 matrix connecting egdes to local edges of elements. 
%
%   QUADRULE is a struct, which specifies a 1D Gauss qaudrature that is
%   used to do the integration on each edge:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   The struct SHAP contains the values and gradients of the shape
%   functions computed by the routine SHAP_HP at the quadrature points of
%   QUADRULE.
%
%   Example:
%
%   p = 10;
%   EDofs = 9*ones(size(Mesh.Edges,1),1);
%   CDofs = 36*ones(size(Mesh.Elements,1),1);
%   Elem2Dof = build_DofMap(Mesh,CDofs,EDofs);
%   QuadRule = Duffy(TProd(gauleg(0,1,2*p)));
%   Shap = shap_hp(QuadRule.x,p);
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemDir_hp(Mesh,Eleme2Dof,L,BdFlags,QuadRule,Shap,FHandle);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nGauss = size(QuadRule.w,1);  % Number of 1D quadrature points
  
  for j1 = BdFlags
      
    % Extract Dirichlet edges
  
    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j1);
    
    for j2 = Loc'
        
      % Compute element map
      
      id_s = Mesh.Edges(j2,1);
      id_e = Mesh.Edges(j2,2);
      Q0 = Mesh.Coordinates(id_s,:);
      Q1 = Mesh.Coordinates(id_e,:);
      dS = norm(Q1-Q0);
      x = ones(nGauss,1)*Q0+QuadRule.x*(Q1-Q0);
      
      % Evaluate Neumann boundary data
         
      FVal = FHandle(x,j1,varargin{:});
      
      % Extract dof numbers on the current edge   
         
      if(Mesh.Edge2Elem(j2,1))
        Elem = Mesh.Edge2Elem(j2,1); 
        ELoc = Mesh.EdgeLoc(j2,1);
      else
        Elem = Mesh.Edge2Elem(j2,2);
        ELoc = Mesh.EdgeLoc(j2,2);
      end
      ndofs = Elem2Dof.EDofs{ELoc}.nDofs(Elem);
      dofs = Elem2Dof.EDofs{ELoc}.Dofs{Elem};
      
      % Compute L2 interpolation for interior edge shape functions 
      
      for i = 1:ndofs
        shap = Shap.eshap{3}.shap{i};
        L(dofs(i)) = L(dofs(i)) + sum(QuadRule.w.*FVal.*shap)*dS;
      end
            
    end
  end
  
return