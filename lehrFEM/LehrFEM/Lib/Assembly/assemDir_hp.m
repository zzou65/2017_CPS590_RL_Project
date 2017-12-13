function [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,BdFlags,QuadRule,Shap,FHandle,varargin)
% ASSEMDIR_HP Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_HP(MESH,ELEM2DOF,BDFLAGS,QUADRULE,SHAP,FHANDLE)
%   incoporates the Dirichlet boundary conditions with the data given by
%   FHANDLE into the finite element solution U. 
%
%   The struct ELEM2DOF describes the element to dof mapping. The boundary
%   condition is only enforced at the vertices of the edges whose boundary
%   flag is equal to the integer BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_HP(MESH,ELEM2DOF,BDFLAGS,QUADRULE,SHAP,FHANDLE, ... 
%   FPARAM) also handles the variable length argument list FPARAM to the
%   boundary data function FHANDLE.
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
%    EDGELOC     P-by-2 matrix connecting egdes to local edges of elements. 
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
%   FREEDOFS is a Q-by-1 matrix specifying the dofs with no prescribed
%   Dirichlet boundary data.  
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
%   [U,FreeDofs] = assemDir_hp(Mesh,Elem2Dof,BdFlags,QuadRule,Shap,FHandle);

%   Copyright 2006-2009 Patrick Meury & Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nGauss = size(QuadRule.w,1);
  
  tmp = [];
  U = zeros(nCoordinates+Elem2Dof.tot_EDofs+Elem2Dof.tot_CDofs,1); 
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
      x = ones(nGauss,1)*Q0+QuadRule.x*(Q1-Q0);
      
      
      % Evaluate Dirichlet boundary data
      
      U(id_s) = FHandle(Q0,j1,varargin{:});
      U(id_e) = FHandle(Q1,j1,varargin{:}); 
      FVal = FHandle(x,j1,varargin{:}) - ...
             U(id_s)*Shap.vshap.shap{1} - ...
             U(id_e)*Shap.vshap.shap{2};
      
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
      Dir = -1; % wrong computation of the direction flag throughout the code
      
      % Compute L2 interpolation for interior edge shape functions 
      
      M = zeros(ndofs);
      L = zeros(ndofs,1);
      for i = 1:ndofs
        shap_1 = Dir^(i+1)*Shap.eshap{3}.shap{i};
        L(i) = sum(QuadRule.w.*FVal.*shap_1);
        for j = 1:ndofs
          shap_2 = Dir^(j+1)*Shap.eshap{3}.shap{j};
          M(i,j) = sum(QuadRule.w.*shap_1.*shap_2);
        end
      end
      tri = triu(M);
      M = tril(transpose(tri),-1)+tri;
      U(dofs) = M\L;
      
      % Collect Dirichlet dofs in temporary container
      
      tmp = [tmp; id_s; id_e; dofs'];
      
    end
  end
  
  % Compute set of free dofs
  
  dofs = 1:(nCoordinates+Elem2Dof.tot_EDofs+Elem2Dof.tot_CDofs);
  FreeDofs = setdiff(dofs,unique(tmp));
  
return