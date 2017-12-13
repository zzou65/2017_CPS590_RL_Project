function L = assemLoad_hp(Mesh,Elem2Dof,QuadRule,Shap,FHandle,varargin)
% ASSEMLOAD_HP Assemble hp-FEM contributions.
%
%   L = ASSEMLOAD_HP(MESH,ELEM2DOF,QUADRULE,SHAP,FHANDLE) assembles the
%   global load vector for the load data given by the function handle
%   FHANDLE. The struct ELEM2DOF describes the element to dof mapping. 
%
%   The struct MESH must at least contain the following fields:
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
%   QUADRULE is a struct, which specifies a 2D Gauss qaudrature that is
%   used to do the integration on each element:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   The struct SHAP contains the values and gradients of the shape
%   functions computed by the routine SHAP_HP at the quadrature points of
%   QUADRULE.
%
%   L = ASSEMLOAD_HP(MESH,ELEM2DOF,QUADRULE,SHAP,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
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
%   L = assemLoad_LFE(Mesh,Elem2Dof,QuadRule,Shap,FHandle);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);  % Number of vertices
  nElements = size(Mesh.Elements,1);        % Number of element of the mesh
  
  % Assemble load vector
  
  L = zeros(nCoordinates+Elem2Dof.tot_EDofs+Elem2Dof.tot_CDofs,1);
  
  % Assemble global load vector
  
  EDofs = zeros(1,3);
  EDir = zeros(1,3);
  for i = 1:nElements
  
    % Extract vertices of the current element
    
    vidx = Mesh.Elements(i,:);
    
    % Extract local polynomial orders
    
    EDofs(1) = Elem2Dof.EDofs{1}.nDofs(i);
    EDofs(2) = Elem2Dof.EDofs{2}.nDofs(i);
    EDofs(3) = Elem2Dof.EDofs{3}.nDofs(i);
    CDofs = Elem2Dof.CDofs.nDofs(i);
    
    % Extract local edge orientations
    
    EDir(1) = Elem2Dof.EDofs{1}.Dir(i);
    EDir(2) = Elem2Dof.EDofs{2}.Dir(i);
    EDir(3) = Elem2Dof.EDofs{3}.Dir(i);
      
    % Preallocate memory for element load vector
  
    Lloc = zeros(3+sum(EDofs)+CDofs,1);
  
    % Compute element map
  
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
  
    x = QuadRule.x*BK+ones(size(QuadRule.w,1),1)*bK;

    % Compute function values
  
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
  
    % Compute element load vector for vertex shape functions
  
    for j = 1:3
      shap = Shap.vshap.shap{j};
      Lloc(j) = sum(QuadRule.w.*FVal.*shap)*det_BK;
    end
  
    % Compute element load vector for edge shape functions
  
    offset = 3;
    for j = 1:EDofs(1)
      shap = EDir(1)^(j+1)*Shap.eshap{1}.shap{j};
      Lloc(j+offset) = sum(QuadRule.w.*FVal.*shap)*det_BK;    
    end
  
    offset = 3+EDofs(1);
    for j = 1:EDofs(2)
      shap = EDir(2)^(j+1)*Shap.eshap{2}.shap{j};
      Lloc(j+offset) = sum(QuadRule.w.*FVal.*shap)*det_BK;  
    end
  
    offset = 3+EDofs(1)+EDofs(2);
    for j = 1:EDofs(3)
      shap = EDir(3)^(j+1)*Shap.eshap{3}.shap{j};
      Lloc(j+offset) = sum(QuadRule.w.*FVal.*shap)*det_BK;  
    end
  
    % Compute element load vector for element shape functions
  
    offset = 3+sum(EDofs);
    for j = 1:CDofs
      shap = Shap.cshap.shap{j};
      Lloc(j+offset) = sum(QuadRule.w.*FVal.*shap)*det_BK;  
    end
    
    % Add local contributions of vertex dofs to global load vector
    
    L(vidx) = L(vidx) + Lloc(1:3);
         
    % Add local contributions of edge dofs to global load vector
    
    loc = 3+(1:EDofs(1));
    L(Elem2Dof.EDofs{1}.Dofs{i}) = L(Elem2Dof.EDofs{1}.Dofs{i})+Lloc(loc);  
      
    loc = 3+EDofs(1)+(1:EDofs(2));    
    L(Elem2Dof.EDofs{2}.Dofs{i}) = L(Elem2Dof.EDofs{2}.Dofs{i})+Lloc(loc);  
      
    loc = 3+EDofs(1)+EDofs(2)+(1:EDofs(3));
    L(Elem2Dof.EDofs{3}.Dofs{i}) = L(Elem2Dof.EDofs{3}.Dofs{i})+Lloc(loc);  
               
    % Add local contributions of element dofs to global load vector
    
    loc = 3+sum(EDofs)+(1:CDofs);
    L(Elem2Dof.CDofs.Dofs{i}) = L(Elem2Dof.CDofs.Dofs{i})+Lloc(loc);  
    
  end
  
return
