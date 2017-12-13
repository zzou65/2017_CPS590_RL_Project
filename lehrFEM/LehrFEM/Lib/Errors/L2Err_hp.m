function err = L2Err_hp(Mesh,u,Elem2Dof,QuadRule,Shap,FHandle,varargin)
% L2ERR_HP Discretization error in L2 norm for hp finite elements.
%
%   ERR = L2ERR_HP(MESH,U,ELEM2DOF,QUADRULE,SHAP,FHANDLE) computes the
%   discretization error between the exact solution given by the function
%   handle FHANDLE and the finite element solution U on the struct MESH.
%
%   The struct ELEM2DOF specifies the element to dof mapping. 
%
%   The struct MESH must at least contain the following fields:
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
%   The struct SHAP contains the values and gradients of the shape
%   functions computed by the routine SHAP_HP at the quadrature points of
%   QUADRULE.
%
%   ERR = L2ERR_HP(MESH,U,ELEM2DOF,QUADRULE,SHAP,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the exact solution
%   FHANDLE.
%
%   Example:
% 
%   err = L2Err_hp(Mesh,u,Elem2Dof,QuadRule,Shap,Fhandle);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nPts = size(QuadRule.w,1);
  
  err = 0;
  for i = 1:nElements
      
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
      
    % Compute element map
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    % Transform quadrature points
    
    x = QuadRule.x*BK+ones(nPts,1)*bK;;
    
    % Evaluate exact and finite element solution
    
    u_EX = FHandle(x,varargin{:});
    
    u_FE = zeros(size(u_EX)); 
    for j = 1:3
      u_FE = u_FE + u(vidx(j))*Shap.vshap.shap{j};
    end
    
    ndofs = Elem2Dof.EDofs{1}.nDofs(i);
    dofs = Elem2Dof.EDofs{1}.Dofs{i};
    dir = Elem2Dof.EDofs{1}.Dir(i);
    for j = 1:ndofs
      u_FE = u_FE + u(dofs(j))*dir^(j+1)*Shap.eshap{1}.shap{j};    
    end
    
    ndofs = Elem2Dof.EDofs{2}.nDofs(i);
    dofs = Elem2Dof.EDofs{2}.Dofs{i};
    dir = Elem2Dof.EDofs{2}.Dir(i);
    for j = 1:ndofs
      u_FE = u_FE + u(dofs(j))*dir^(j+1)*Shap.eshap{2}.shap{j};    
    end
    
    ndofs = Elem2Dof.EDofs{3}.nDofs(i);
    dofs = Elem2Dof.EDofs{3}.Dofs{i};
    dir = Elem2Dof.EDofs{3}.Dir(i);
    for j = 1:ndofs
      u_FE = u_FE + u(dofs(j))*dir^(j+1)*Shap.eshap{3}.shap{j};    
    end
    
    ndofs = Elem2Dof.CDofs.nDofs(i);
    dofs = Elem2Dof.CDofs.Dofs{i};
    for j = 1:ndofs
      u_FE = u_FE + u(dofs(j))*Shap.cshap.shap{j};  
    end
    
    % Compute discretization error
    
    err = err + sum(QuadRule.w.*(u_FE-u_EX).^2)*det_BK;
      
  end
  err = sqrt(err);
  
return