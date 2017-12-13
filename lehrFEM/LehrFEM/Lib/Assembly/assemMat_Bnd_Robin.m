function varargout = assemMat_Bnd_Robin(Mesh,EHandle,prec,points)
% ASSEMMAT_BND_Robin Assemble boundary part of laplace with robin boundary 
% condition in LFE.
%
%   A = ASSEMMAT_BND_Robin(MESH,EHANDLE) Assemble boundary part of laplace 
%   with robin boundary condition in LFE given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_BND_Robin(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_BND_Robin(MESH,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGES        N-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-3 matrix connecting egdes to local edges of elements. 
%    
%   Example:
%
%   A = assemMat_Bnd_Robin(Mesh,@STIMA_Bnd_Robin);
%  
%   See also set_Rows, set_Cols.
%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements in the mesh
  nEdges = size(Mesh.Edges,1);        % Number of edges in the mesh
  nCoord = size(Mesh.Coordinates,1);        % Number of edges in the mesh

  QuadRule = gauleg(0,1,points,prec);
  % Allocate memory

  I = zeros(4*nEdges,1);
  J = zeros(4*nEdges,1);
  A = zeros(4*nEdges,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Assemble element contributions
  
  loc = 1:4;
  last=0;
  for i = 1:nEdges
  
    % Check for boundary edge  
      
    if(BdFlags(i) < 0)  
    
        
      vid=Mesh.Edges(i,:);
      Coordinates = Mesh.Coordinates(vid,:);
      Edgevol=norm(Coordinates(1,:)-Coordinates(2,:));
      
      t=QuadRule.x;
      x=t*Coordinates(1,:)+(1-t)*Coordinates(2,:);
      g=EHandle(x);
      shap1=t;
      shap2=1-t;
      % Compute element contributions
      Aloc(1,1)=Edgevol*sum(QuadRule.w.*g.*shap1.*shap1);
      Aloc(1,2)=Edgevol*sum(QuadRule.w.*g.*shap1.*shap2);
      Aloc(2,1)=Aloc(1,2);
      Aloc(2,2)=Edgevol*sum(QuadRule.w.*g.*shap2.*shap2);
    
      % Add element contributions to stiffness matrix
    
      I(loc) = set_Rows(vid,2);
      J(loc) = set_Cols(vid,2);
      A(loc) = Aloc(:);
      loc = loc + 4;
      last=last+4;
     end
      
  end
  
  % Assign output arguments
  loc=1:last;
  if(nargout > 1)
    varargout{1} = I(loc);
    varargout{2} = J(loc);
    varargout{3} = A(loc);
  else
    varargout{1} = sparse(I(loc),J(loc),A(loc),nCoord,nCoord);  
  end
  
return
