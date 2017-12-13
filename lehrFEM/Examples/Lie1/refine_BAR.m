function  New_Mesh = refine_BAR(Old_Mesh,varargin)
% REFINE_BAR Barycentric refinement.
%
%   MESH = REFINE_Bar(MESH) performs one barycentric refinement step with the
%   struct MESH.
%
%   MESH = REFINE_BAR(MESH,DHANDLE) performs onebarycentric refinement step
%   with the struct MESH. During red refinement the signed distance function 
%   DHANDLE (function handle/inline onject) is used to project the new vertices
%   on the boundary edges onto the boundary of the domain.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two vertices
%                 i and j are connected by an edge with number VERT2EDGE(i,j).
%
%   MESH = REFINE_Bar(MESH,DHANDLE,DPARAM) also handles the variable argument
%   list DPARAM to the signed distance function DHANDLE.
%
%   Example:
%
%   Mesh = refine_Bar(Mesh,@dist_circ,[0 0],1);
%
%   See also ADD_EDGES.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nCoordinates = size(Old_Mesh.Coordinates,1);
  [nElements,nVert] = size(Old_Mesh.Elements);
  nEdges = size(Old_Mesh.Edges,1);
  nBdEdges = size(find(Old_Mesh.BdFlags < 0),1);
  
  % Extract input arguments
  
  if(nargin > 1)
    DHandle = varargin{1};
    DParam = varargin(2:nargin-1);
  else
    DHandle = [];
  end
  
  % Red refinement for triangular or triangular elements
  
    % Preallocate memory
        
    New_Mesh.Coordinates = zeros(nCoordinates+nEdges+nElements,2);
    New_Mesh.Elements = zeros(6*nElements,3); 
    New_Mesh.BdFlags = zeros(2*nEdges+6*nElements,1);
    if(isfield(Old_Mesh,'ElemFlag'))
      New_Mesh.ElemFlag = zeros(6*nElements,1);    
    end
    
    % Do barycentric refinement
  
    New_Mesh.Coordinates(1:nCoordinates,:) = Old_Mesh.Coordinates;
    Bd_Idx = 0;
    Aux = zeros(nBdEdges,4);
    for i = 1:nElements
    
      % Get vertex numbers of the current element
      
      i1 = Old_Mesh.Elements(i,1);
      i2 = Old_Mesh.Elements(i,2);
      i3 = Old_Mesh.Elements(i,3);
    
      % Compute vertex numbers of new vertices localized on edges
    
      j1 = nCoordinates+Old_Mesh.Vert2Edge(i2,i3);
      j2 = nCoordinates+Old_Mesh.Vert2Edge(i3,i1);
      j3 = nCoordinates+Old_Mesh.Vert2Edge(i1,i2);
      b=nCoordinates+nEdges+i;
      % Generate new elements
    
      New_Mesh.Elements(6*(i-1)+1,:) = [i1 j3 b];
      New_Mesh.Elements(6*(i-1)+2,:) = [j3 i2 b];
      New_Mesh.Elements(6*(i-1)+3,:) = [i2 j1 b];
      New_Mesh.Elements(6*(i-1)+4,:) = [j1 i3 b];
      New_Mesh.Elements(6*(i-1)+5,:) = [i3 j2 b];
      New_Mesh.Elements(6*(i-1)+6,:) = [j2 i1 b];
      
      % Generate new vertex on edge 1, project to boundary if necessary
     
      BdFlag_1 = Old_Mesh.BdFlags(Old_Mesh.Vert2Edge(i2,i3));
      if(BdFlag_1 < 0)
        if(~isempty(DHandle))
          DEPS = sqrt(eps)*norm(Old_Mesh.Coordinates(i2,:)-Old_Mesh.Coordinates(i3,:));
          x = 1/2*(Old_Mesh.Coordinates(i2,:)+Old_Mesh.Coordinates(i3,:));
          dist = feval(DHandle,x,DParam{:});
          grad_dist = ([feval(DHandle,x+[DEPS 0],DParam{:}) feval(DHandle,x+[0 DEPS],DParam{:})]-dist)/DEPS;  
          New_Mesh.Coordinates(j1,:) = x-dist*grad_dist;
        else
          New_Mesh.Coordinates(j1,:) = 1/2*(Old_Mesh.Coordinates(i2,:)+Old_Mesh.Coordinates(i3,:));  
        end
        Bd_Idx = Bd_Idx+1;
        Aux(Bd_Idx,:) = [BdFlag_1 i2 j1 i3];
      else
        New_Mesh.Coordinates(j1,:) = 1/2*(Old_Mesh.Coordinates(i2,:)+Old_Mesh.Coordinates(i3,:)); 
      end
          
      % Generate new vertex on edge 2, project to boundary if necessary
    
      BdFlag_2 = Old_Mesh.BdFlags(Old_Mesh.Vert2Edge(i3,i1));
      if(BdFlag_2 < 0)
        if(~isempty(DHandle))
          DEPS = sqrt(eps)*norm(Old_Mesh.Coordinates(i3,:)-Old_Mesh.Coordinates(i1,:));  
          x = 1/2*(Old_Mesh.Coordinates(i3,:)+Old_Mesh.Coordinates(i1,:));
          dist = feval(DHandle,x,DParam{:});
          grad_dist = ([feval(DHandle,x+[DEPS 0],DParam{:}) feval(DHandle,x+[0 DEPS],DParam{:})]-dist)/DEPS;  
          New_Mesh.Coordinates(j2,:) = x-dist*grad_dist;
        else
          New_Mesh.Coordinates(j2,:) = 1/2*(Old_Mesh.Coordinates(i3,:)+Old_Mesh.Coordinates(i1,:));   
        end
        Bd_Idx = Bd_Idx+1;
        Aux(Bd_Idx,:) = [BdFlag_2 i3 j2 i1];
      else
        New_Mesh.Coordinates(j2,:) = 1/2*(Old_Mesh.Coordinates(i3,:)+Old_Mesh.Coordinates(i1,:));
      end
    
      % Generate new vertex on edge 3, project to boundary if necessary
    
      BdFlag_3 = Old_Mesh.BdFlags(Old_Mesh.Vert2Edge(i1,i2));
      if(BdFlag_3 < 0)
        if(~isempty(DHandle))
          DEPS = sqrt(eps)*norm(Old_Mesh.Coordinates(i1,:)-Old_Mesh.Coordinates(i2,:));
          x = 1/2*(Old_Mesh.Coordinates(i1,:)+Old_Mesh.Coordinates(i2,:));
          dist = feval(DHandle,x,DParam{:});
          grad_dist = ([feval(DHandle,x+[DEPS 0],DParam{:}) feval(DHandle,x+[0 DEPS],DParam{:})]-dist)/DEPS;  
          New_Mesh.Coordinates(j3,:) = x-dist*grad_dist;
        else
          New_Mesh.Coordinates(j3,:) = 1/2*(Old_Mesh.Coordinates(i1,:)+Old_Mesh.Coordinates(i2,:));  
        end
        Bd_Idx = Bd_Idx+1;
        Aux(Bd_Idx,:) = [BdFlag_3 i1 j3 i2];
      else
        New_Mesh.Coordinates(j3,:) = 1/2*(Old_Mesh.Coordinates(i1,:)+Old_Mesh.Coordinates(i2,:));
      end
      
      % Generate new vertex at barycenter of element i
      New_Mesh.Coordinates(b,:) = 1/3*(Old_Mesh.Coordinates(i1,:)+Old_Mesh.Coordinates(i2,:)+Old_Mesh.Coordinates(i3,:)); 
      % Update element flag if defined
      
      if(isfield(Old_Mesh,'ElemFlag'))
        New_Mesh.ElemFlag(6*(i-1)+1) = Old_Mesh.ElemFlag(i);
        New_Mesh.ElemFlag(6*(i-1)+2) = Old_Mesh.ElemFlag(i);
        New_Mesh.ElemFlag(6*(i-1)+3) = Old_Mesh.ElemFlag(i);
        New_Mesh.ElemFlag(6*(i-1)+4) = Old_Mesh.ElemFlag(i);
        New_Mesh.ElemFlag(6*(i-1)+5) = Old_Mesh.ElemFlag(i);
        New_Mesh.ElemFlag(6*(i-1)+6) = Old_Mesh.ElemFlag(i);
      end
    end
  
    % Add edges to new mesh
    New_Mesh.Max_Nodes=2*Old_Mesh.Max_Nodes;
    New_Mesh = add_Edges(New_Mesh);
  
    % Update boundary flags
  
    for i = 1:nBdEdges
      New_Mesh.BdFlags(New_Mesh.Vert2Edge(Aux(i,2),Aux(i,3))) = Aux(i,1);
      New_Mesh.BdFlags(New_Mesh.Vert2Edge(Aux(i,3),Aux(i,4))) = Aux(i,1);
    end
    
return
