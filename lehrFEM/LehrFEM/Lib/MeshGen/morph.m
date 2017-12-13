function New_Mesh = morph(Old_Mesh)
% MORPH Morph triangular into quadrilateral meshes.
%
%   MESH = MORPH(MESH) morphs any 2D triangular mesh into a 2D quadrilateral
%   mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   Mesh = morph(Mesh);
%
%   See also ADD_EDGES, ADD_EDGE2ELEM, GET_BDEDGES, QUADQUAL.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  QMIN = 0.5;

  nCoordinates = size(Old_Mesh.Coordinates,1);
  nElements = size(Old_Mesh.Elements,1);
  
  % Add edge data structure to mesh

  Old_Mesh = add_Edges(Old_Mesh);
  Old_Mesh = add_Edge2Elem(Old_Mesh);
  BdLoc = get_BdEdges(Old_Mesh);

  % Mark all non-removable edges
  
  nEdges = size(Old_Mesh.Edges,1);
  isdiagonal = ones(nEdges,1);
  isdiagonal(BdLoc) = 0;
    
  % Compute distortion factors and precompute quadrilaterals
  
  beta = zeros(nEdges,1);
  Quad = zeros(nEdges,4);
  for i = 1:nEdges
    if(isdiagonal(i))
       ElemLeft = Old_Mesh.Edge2Elem(i,1);
       ElemRight = Old_Mesh.Edge2Elem(i,2);
       OppLeft = Old_Mesh.EdgeLoc(i,1);
       OppRight = Old_Mesh.EdgeLoc(i,2);
       Quad(i,:) = [Old_Mesh.Elements(ElemRight,OppRight), ...
                    Old_Mesh.Elements(ElemRight,rem(OppRight+3,3)+1), ...
                    Old_Mesh.Elements(ElemLeft,OppLeft), ...
                    Old_Mesh.Elements(ElemLeft,rem(OppLeft+3,3)+1)];
       Coord = Old_Mesh.Coordinates(Quad(i,:),:);
       beta(i) = QuadQual(Coord);    
    end
  end
  
  % Sort distortion factors 
  
  [beta_s,loc] = sort(beta);
  
  % Merge triangles into quadrilaterals
  
  tmp = zeros(3*nElements,4);
  isdeleted = zeros(nEdges,1);
  ismerged = zeros(nElements,1);
  nQuad = 0;
  i = 0;
  for j = size(beta_s,1):-1:1
    if(isdiagonal(loc(j)) && beta_s(j) > QMIN)
      i = loc(j);
      break
    else  
      beta(j) = [];
      loc(j) = [];
    end
  end
  while(i)
    nQuad = nQuad+1;  
    
    % Remove edge and merge triangles
    
    isdiagonal(i) = 0;
    isdeleted(i) = 1;
    ismerged(Old_Mesh.Edge2Elem(i,1)) = 1;
    ismerged(Old_Mesh.Edge2Elem(i,2)) = 1;
    tmp(nQuad,:) = Quad(i,:);
    
    % Update diagonal flags
    
    isdiagonal(Old_Mesh.Vert2Edge(Quad(i,1),Quad(i,2))) = 0;
    isdiagonal(Old_Mesh.Vert2Edge(Quad(i,2),Quad(i,3))) = 0;
    isdiagonal(Old_Mesh.Vert2Edge(Quad(i,3),Quad(i,4))) = 0;
    isdiagonal(Old_Mesh.Vert2Edge(Quad(i,4),Quad(i,1))) = 0;   
    
    % Get next removable diagonal
    
    i = 0;
    for j = size(beta_s,1):-1:1
      if(isdiagonal(loc(j)) && beta_s(j) > QMIN)
        i = loc(j);
        break
      else
        beta_s(j) = [];
        loc(j) = [];
      end
    end
  end
  Quad = tmp(1:nQuad,:);
  tmp = 0;
  
  % Allocate memory for new mesh
  
  New_Mesh.Coordinates = zeros(nCoordinates+nEdges+nElements,2);
  New_Mesh.Elements = zeros(3*nElements,4);   
  nCoord = nCoordinates;
  New_Mesh.Coordinates(1:nCoord,:) = Old_Mesh.Coordinates; 
  
  % Create new vertices on edges
  
  cnt = 0;
  VertOnEdge = zeros(nEdges,1);
  for i = 1:nEdges
    if(~isdeleted(i))
      cnt = cnt+1;
      VertOnEdge(i) = nCoord+cnt;
      New_Mesh.Coordinates(nCoord+cnt,:) = 1/2*(New_Mesh.Coordinates(Old_Mesh.Edges(i,1),:) + ...
                                                New_Mesh.Coordinates(Old_Mesh.Edges(i,2),:));                                  
    end
  end
  nCoord = nCoord+cnt;
          
  % Create quadrilateral elements  
    
  nElem = 0;
  
  % Create new quadrilaterals from merged triangles
  
  for i = 1:nQuad
      
    % Create new vertices inside quadrilaterals  
    
    jc = nCoord+i;
    New_Mesh.Coordinates(jc,:) = 1/4*(New_Mesh.Coordinates(Quad(i,1),:) + ...
                                      New_Mesh.Coordinates(Quad(i,2),:) + ...
                                      New_Mesh.Coordinates(Quad(i,3),:) + ...
                                      New_Mesh.Coordinates(Quad(i,4),:));
    
    % Get vertices of current merged triangles
    
    i1 = Quad(i,1);
    i2 = Quad(i,2);
    i3 = Quad(i,3);
    i4 = Quad(i,4);
    
    j1 = VertOnEdge(Old_Mesh.Vert2Edge(i1,i2));
    j2 = VertOnEdge(Old_Mesh.Vert2Edge(i2,i3));
    j3 = VertOnEdge(Old_Mesh.Vert2Edge(i3,i4));
    j4 = VertOnEdge(Old_Mesh.Vert2Edge(i4,i1));
                                  
    % Create new quadrilaterals
                                        
    nElem = nElem+1;
    New_Mesh.Elements(nElem,:) = [i1 j1 jc j4];
    nElem = nElem+1;
    New_Mesh.Elements(nElem,:) = [j1 i2 j2 jc];
    nElem = nElem+1;
    New_Mesh.Elements(nElem,:) = [j4 jc j3 i4];
    nElem = nElem+1;
    New_Mesh.Elements(nElem,:) = [jc j2 i3 j3];
  end
  nCoord = nCoord+nQuad;
  
  % Create new quadrilaterals from leftover tringles
  
  cnt = 0;
  for i = 1:nElements
    if(~ismerged(i))
      cnt = cnt+1;
        
      % Create new vertices inside triangles
      
      jc = nCoord+cnt;
      New_Mesh.Coordinates(jc,:) = 1/3*(New_Mesh.Coordinates(Old_Mesh.Elements(i,1),:) + ...
                                        New_Mesh.Coordinates(Old_Mesh.Elements(i,2),:) + ...
                                        New_Mesh.Coordinates(Old_Mesh.Elements(i,3),:));
      
      % Get vertices of current triangle                              
         
      i1 = Old_Mesh.Elements(i,1);
      i2 = Old_Mesh.Elements(i,2);
      i3 = Old_Mesh.Elements(i,3);
      
      j1 = VertOnEdge(Old_Mesh.Vert2Edge(i2,i3));
      j2 = VertOnEdge(Old_Mesh.Vert2Edge(i3,i1));
      j3 = VertOnEdge(Old_Mesh.Vert2Edge(i1,i2));
      
      % Create new quadrilaterals                                      
                                            
      nElem = nElem+1;
      New_Mesh.Elements(nElem,:) = [i1 j3 jc j2];
      nElem = nElem+1;
      New_Mesh.Elements(nElem,:) = [j3 i2 j1 jc];
      nElem = nElem+1;
      New_Mesh.Elements(nElem,:) = [j2 jc j1 i3];
    end
  end
  nCoord = nCoord+cnt;
  
  % Squeeze element and vertex lists
  
  New_Mesh.Coordinates = New_Mesh.Coordinates(1:nCoord,:);
  New_Mesh.Elements = New_Mesh.Elements(1:nElem,:);
  
return
  