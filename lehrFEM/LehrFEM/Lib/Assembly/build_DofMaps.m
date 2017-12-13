function Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs)
% BUILD_DOFMAPS Generate Element to Dof mapping.
%
%   ELEM2DOF = BUILD_DOFMAPS(MESH,EDOFS,CDOFS) generates the element to dof
%   mapping for the struct mesh and the element and edge polynomial degrees
%   specified by the N-by-1 and P-by-1 matrices CDOFS and EDOFS.
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
%   Example:
%
%   Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);  % Number of vertices
  nElements = size(Mesh.Elements,1);        % Number of elements
  nEdges = size(Mesh.Edges,1);              % Number of edges
  nEDofs = sum(EDofs);                      % Total number of edge dofs      
  nCDofs = sum(CDofs);                      % Total number of element dofs
  
  % Preallocate memory
  
  for i = 1:3
    Elem2Dof.EDofs{i}.Dofs = cell(nElements,1);
    Elem2Dof.EDofs{i}.nDofs = zeros(nElements,1); 
    Elem2Dof.EDofs{i}.Dir = zeros(nElements,1);
  end
  Elem2Dof.CDofs.Dofs = cell(nElements,1);
  Elem2Dof.CDofs.nDofs = zeros(nElements,1);
  Elem2Dof.tot_EDofs = nEDofs;
  Elem2Dof.tot_CDofs = nCDofs;
  
  % Build edge dof map
  
  offset_e = nCoordinates;
  for i = 1:nEdges   
    
    % Build edge dof maps for left hand side element
    
    j = Mesh.Edge2Elem(i,1);
    if(j > 0)
      loc = Mesh.EdgeLoc(i,1);
      Elem2Dof.EDofs{loc}.Dofs{j} = offset_e+(1:EDofs(i));
      Elem2Dof.EDofs{loc}.nDofs(j) = EDofs(i);
      id_s = Mesh.Elements(j,rem(loc+4,3)+1);
      if(id_s == Mesh.Edges(i,1))
        Elem2Dof.EDofs{loc}.Dir(j) = 1;
      else
        Elem2Dof.EDofs{loc}.Dir(j) = -1;
      end
    end
    
    % Build dof maps for right hand side element
    
    j = Mesh.Edge2Elem(i,2);
    if(j > 0) 
      loc = Mesh.EdgeLoc(i,2); 
      Elem2Dof.EDofs{loc}.Dofs{j} = offset_e+(1:EDofs(i));
      Elem2Dof.EDofs{loc}.nDofs(j) = EDofs(i);
      id_s = Mesh.Elements(j,rem(loc+4,3)+1);
      if(id_s == Mesh.Edges(i,1))
        Elem2Dof.EDofs{loc}.Dir(j) = 1; 
      else
        Elem2Dof.EDofs{loc}.Dir(j) = -1; 
      end
    end
           
    offset_e = offset_e+EDofs(i);
  end
  
  % Build element dof map
  
  offset_c = nCoordinates+nEDofs;
  for i = 1:nElements
     Elem2Dof.CDofs.Dofs{i} = offset_c+(1:CDofs(i));
     offset_c = offset_c+CDofs(i);
  end
  Elem2Dof.CDofs.nDofs = CDofs;

return