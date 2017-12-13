function varargout = assemMat_hp(Mesh,Elem2Dof,EHandle,varargin)
% ASSEMMAT_HP Assemble hp-FEM contributions.
%
%   A = ASSEMMAT_HP(MESH,ELEM2DOF,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in a sparse representation. The struct ELEM2DOF 
%   describes the element to dof mapping obtained from the routine
%   BUILD_DOFMAPS.
%
%   A = ASSEMMAT_HP(MESH,ELEM2DOF,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process. 
%
%   [I,J,A] = ASSEMMAT_HP(MESH,ELEM2DOF,EHANDLE) assembles the global
%   matrix from the local element contributions given by the function
%   handle EHANDLE and returns the matrix in an array representation.
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
%   Example:
%
%   p = 10;
%   EDofs = 9*ones(size(Mesh.Edges,1),1);
%   CDofs = 36*ones(size(Mesh.Elements,1),1);
%   Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);
%   QuadRule = Duffy(TProd(gauleg(0,1,2*p)));
%   Shap = shap_hp(QuadRule.x,p);
%   A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule,Shap);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements
                      
  % Preallocate memory
    
  nEntries = sum((3*ones(nElements,1) + Elem2Dof.EDofs{1}.nDofs + ...
                  Elem2Dof.EDofs{2}.nDofs + Elem2Dof.EDofs{3}.nDofs + ...
                  Elem2Dof.CDofs.nDofs).^2);
  I = zeros(nEntries,1);
  J = zeros(nEntries,1);
  A = zeros(nEntries,1);
    
  % Assemble global stiffness matrix  
    
  EDofs = zeros(1,3);
  EDir = zeros(1,3);
  CDofs = 0;
  offset = 0;
  for i = 1:nElements
  
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Extract local polynomial orders
    
    EDofs(1) = Elem2Dof.EDofs{1}.nDofs(i);
    EDofs(2) = Elem2Dof.EDofs{2}.nDofs(i);
    EDofs(3) = Elem2Dof.EDofs{3}.nDofs(i);
    CDofs = Elem2Dof.CDofs.nDofs(i);
      
    % Extract local edge orientations
    
    EDir(1) = Elem2Dof.EDofs{1}.Dir(i);
    EDir(2) = Elem2Dof.EDofs{2}.Dir(i);
    EDir(3) = Elem2Dof.EDofs{3}.Dir(i);
    
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),EDofs,EDir,CDofs,varargin{:});
    
    % Add contributions to global matrix
    
    idx = [vidx ...
           Elem2Dof.EDofs{1}.Dofs{i} ...
           Elem2Dof.EDofs{2}.Dofs{i} ...
           Elem2Dof.EDofs{3}.Dofs{i} ...
           Elem2Dof.CDofs.Dofs{i}];
    n_idx = 3+sum(EDofs)+CDofs;
    
    loc = offset + (1:n_idx^2); 
    I(loc) = set_Rows(idx,n_idx);
    J(loc) = set_Cols(idx,n_idx);
    A(loc) = Aloc(:);
    
    offset = offset + n_idx^2;
    
  end  
      
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);    
  end
  
return