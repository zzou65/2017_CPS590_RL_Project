function [Mesh,fatherVert] = refine_LEB(Mesh,Marked_Elements)
% REFINE_LEB largest edge bisection.
%
%   MESH = REFINE_LEB(MESH,MARKED_ELEMENTS) adaptively refines the struct
%   MESH by using the largest edge bisection on the triangular elements
%   MARKED_ELEMENTS.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES N-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh. 
%    NEIGH       M-by-3 matrix specifying all elements sharing a common
%                edge with the current element. If an edge is part of the
%                boundary then NEIGH contains the boundary flag of that
%                edge.
%    OPP_VERT    M-by-3 matrix specifying the local number of the opposite
%                vertex of the element sharing an edge with the current
%                element.
%
%   [MESH,FATHERVERT] = REFINE_LEB(...) also returns an nx2 matrix
%   FATHERVERT, where n is the number of new vertices.  If the original
%   mesh has N vertices, then FATHERVERT(i,:) stores the endpoints of the
%   edge on which the N+i-th vertex of the refined mesh was constructed.
%   These will always be vertices of the original mesh.
%
%   Example:
%
%   Mesh = init_LEB(Mesh);
%   Mesh = refine_LEB(Mesh,[1 2 4 5 10]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Recursively refine the new mesh
  
  fatherVert = zeros(0,2);
  
  while(~isempty(Marked_Elements))
    Elem = Marked_Elements(1);
    [Mesh,Marked_Elements,fatherVert] = recursive_refine(Elem,Mesh,Marked_Elements,fatherVert);
  end

return


%%% Recursive refinement %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mesh,Marked_Elements,fatherVert] = recursive_refine(Elem,Mesh,Marked_Elements,fatherVert)
% RECURSIVE_REFINE Recursive green refinements.
%
%   [MESH,MARKED_ELEMENTS,FATHERVERT] = RECURSIVE_REFINE(ELEM,MESH,MARKED_ELEMENTS,FATHERVERT)
%   recursively refines the element ELEM from the struct MESH by green
%   refinements.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES N-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh. 
%    NEIGH       M-by-3 matrix specifying all elements sharing a common
%                edge with the current element. If an edge is part of the
%                boundary then NEIGH contains the boundary flag of that
%                edge.
%    OPP_VERT    M-by-3 matrix specifying the local number of the opposite
%                vertex of the element sharing an edge with the current
%                element.
%
%   Example:
%
%   [Mesh,Marked_Elements,fatherVert] = recursive_refine(Elem,Mesh,Marked_Elements,zeros(0,2));
%

%   Copyright 2005-2005 Patrick Meury & Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  Neigh = Mesh.Neigh(Elem,3);
  if(Neigh < 0)
     nCoordinates = size(Mesh.Coordinates,1);
     nElements = size(Mesh.Elements,1); 
      
     % Generate new vertex on boundary
     
     New_Vert = nCoordinates+1; 
     Mesh.Coordinates(New_Vert,:) = 1/2*(Mesh.Coordinates(Mesh.Elements(Elem,1),:)+Mesh.Coordinates(Mesh.Elements(Elem,2),:));
     
     % Store parent vertices of new vertex, ie the endpoints of the edge
     
     fatherVert = [fatherVert; Mesh.Elements(Elem,1:2)];
      
     % Generate new elements
     
     Elem_1 = Elem;
     Elem_2 = nElements+1;
     Mesh.Elements(Elem_2,:) = [Mesh.Elements(Elem,2) Mesh.Elements(Elem,3) New_Vert];
     Mesh.Elements(Elem_1,:) = [Mesh.Elements(Elem,3) Mesh.Elements(Elem,1) New_Vert];
     
     % Update element flag if defined
     
     if(isfield(Mesh,'ElemFlag'))
       Mesh.ElemFlag(Elem_2) = Mesh.ElemFlag(Elem_1);   
     end
     
     % Update neighbours and opposite vertices of neighbour elements
     
     if(Mesh.Neigh(Elem,2) > 0)
       Mesh.Neigh(Mesh.Neigh(Elem,2),Mesh.Opp_Vert(Elem,2)) = Elem_1;
       Mesh.Opp_Vert(Mesh.Neigh(Elem,2),Mesh.Opp_Vert(Elem,2)) = 3;
     end
     if(Mesh.Neigh(Elem,1) > 0)
       Mesh.Neigh(Mesh.Neigh(Elem,1),Mesh.Opp_Vert(Elem,1)) = Elem_2;
       Mesh.Opp_Vert(Mesh.Neigh(Elem,1),Mesh.Opp_Vert(Elem,1)) = 3;
     end
     
     % Update neighbours and opposite vertices of new elements
          
     Mesh.Neigh(Elem_2,:) = [Elem_1 Mesh.Neigh(Elem,3) Mesh.Neigh(Elem,1)];
     Mesh.Neigh(Elem_1,:) = [Mesh.Neigh(Elem,3) Elem_2 Mesh.Neigh(Elem,2)]; 
     Mesh.Opp_Vert(Elem_2,:) = [2 0 Mesh.Opp_Vert(Elem,1)];
     Mesh.Opp_Vert(Elem_1,:) = [0 1 Mesh.Opp_Vert(Elem,2)];
     
     % Remove refined elements from list of marked elements
  
     Marked_Elements = Marked_Elements(Marked_Elements ~= Elem);
  else
          
    % Check whether neighbour element is compatibly divisible
      
    if(Mesh.Opp_Vert(Elem,3)~= 3)
      [Mesh,Marked_Elements,fatherVert] = recursive_refine(Neigh,Mesh,Marked_Elements,fatherVert);
    end
    Neigh = Mesh.Neigh(Elem,3);
    nCoordinates = size(Mesh.Coordinates,1);
    nElements = size(Mesh.Elements,1);
    
    % Generate new vertex
    
    New_Vert = nCoordinates+1;
    Mesh.Coordinates(New_Vert,:) = 1/2*(Mesh.Coordinates(Mesh.Elements(Elem,1),:)+Mesh.Coordinates(Mesh.Elements(Elem,2),:));
  
    % Store parent vertices of new vertex, ie the endpoints of the edge
    
    fatherVert = [fatherVert; Mesh.Elements(Elem,1:2)];
    
    % Generate new elements
  
    Elem_1 = Elem;
    Elem_2 = nElements+1;
    Elem_3 = Neigh;
    Elem_4 = nElements+2;
    Mesh.Elements(Elem_2,:) = [Mesh.Elements(Elem,2) Mesh.Elements(Elem,3) New_Vert];
    Mesh.Elements(Elem_1,:) = [Mesh.Elements(Elem,3) Mesh.Elements(Elem,1) New_Vert];
    Mesh.Elements(Elem_4,:) = [Mesh.Elements(Neigh,3) Mesh.Elements(Neigh,1) New_Vert];
    Mesh.Elements(Elem_3,:) = [Mesh.Elements(Neigh,2) Mesh.Elements(Neigh,3) New_Vert];
  
    % Update element flag if defined
    
    if(isfield(Mesh,'ElemFlag'))
      Mesh.ElemFlag(Elem_2) = Mesh.ElemFlag(Elem_1);
      Mesh.ElemFlag(Elem_4) = Mesh.ElemFlag(Elem_3);
    end
    
    % Update neighbours and opposite vertices of neighbour elements
    
    if(Mesh.Neigh(Elem,1) > 0)
      Mesh.Neigh(Mesh.Neigh(Elem,1),Mesh.Opp_Vert(Elem,1)) = Elem_2;
      Mesh.Opp_Vert(Mesh.Neigh(Elem,1),Mesh.Opp_Vert(Elem,1)) = 3;
    end
    if(Mesh.Neigh(Elem,2) > 0)
      Mesh.Neigh(Mesh.Neigh(Elem,2),Mesh.Opp_Vert(Elem,2)) = Elem_1;
      Mesh.Opp_Vert(Mesh.Neigh(Elem,2),Mesh.Opp_Vert(Elem,2)) = 3;
    end
    if(Mesh.Neigh(Neigh,2) > 0)
      Mesh.Neigh(Mesh.Neigh(Neigh,2),Mesh.Opp_Vert(Neigh,2)) = Elem_4;
      Mesh.Opp_Vert(Mesh.Neigh(Neigh,2),Mesh.Opp_Vert(Neigh,2)) = 3;
    end
    if(Mesh.Neigh(Neigh,1) > 0)
      Mesh.Neigh(Mesh.Neigh(Neigh,1),Mesh.Opp_Vert(Neigh,1)) = Elem_3;
      Mesh.Opp_Vert(Mesh.Neigh(Neigh,1),Mesh.Opp_Vert(Neigh,1)) = 3;
    end    
    
    % Update neighbours of new elements
    
    Mesh.Neigh(Elem_2,:) = [Elem_1 Elem_4 Mesh.Neigh(Elem,1)];
    Mesh.Neigh(Elem_1,:) = [Elem_3 Elem_2 Mesh.Neigh(Elem,2)];
    Mesh.Neigh(Elem_4,:) = [Elem_2 Elem_3 Mesh.Neigh(Neigh,2)];
    Mesh.Neigh(Elem_3,:) = [Elem_4 Elem_1 Mesh.Neigh(Neigh,1)];
    Mesh.Opp_Vert(Elem_2,:) = [2 1 Mesh.Opp_Vert(Elem,1)];
    Mesh.Opp_Vert(Elem_1,:) = [2 1 Mesh.Opp_Vert(Elem,2)];
    Mesh.Opp_Vert(Elem_4,:) = [2 1 Mesh.Opp_Vert(Neigh,2)];
    Mesh.Opp_Vert(Elem_3,:) = [2 1 Mesh.Opp_Vert(Neigh,1)];
    
    % Remove element from list of marked elements
  
    Marked_Elements = Marked_Elements(Marked_Elements ~= Elem & Marked_Elements ~= Neigh);
  end  
  
return