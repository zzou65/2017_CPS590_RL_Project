function New_Mesh = add_AdjNodes(Old_Mesh)
% ADD_AdjNodes Adds ids of adjacent nodes to the mesh.
%
%   MESH = ADD_AdjNodes(MESH) adds adjacent nodes to the mesh:
%    ADJNODES: M-by-Q matrix specifying the nodes adjacent to a vertex with 
%              Coordinates(i), in sequential order
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%    ADJELEMENTS  M-by-Q matrix specifying the elements of the mesh sharing
%                 vertex i of the mesh.
%    NADJELEMENTS M-by-1 matrix specifying the actual number of elements
%                 sharing vertex COORDINATES(i) of the mesh.  
%
%   Example:
%
%   Mesh = add_Patches(Mesh);

%   Copyright 2006-2006 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   nCoordinates = size(Old_Mesh.Coordinates,1);
   [nElements,nVert] = size(Old_Mesh.Elements);
   Old_Mesh=add_Patches(Old_Mesh);
   % Initialize constants
   
   MAX_ELEMENTS = max(Old_Mesh.nAdjElements);
   
   % Preallocate memory
   
   AdjNodes = zeros(nCoordinates,MAX_ELEMENTS);
   
   % Build list of adjacent elements
   
   for i = 1:nCoordinates
       nodes=Old_Mesh.Elements(setdiff(Old_Mesh.AdjElements(i,:),0),:);
       % check if node i is on boundary and move the nodes of a boundary
       % element on the firt lines od array nodes
       u_nod=unique(nodes);
       if (size(nodes,1)==1)
           AdjNodes(i,[1 2])=setdiff(nodes,i);
       else
        bdflag=find(sum(histc(nodes(:),u_nod),2)==1,2,'first'); 
        if (~isempty(bdflag))
           [m n]=find(nodes==u_nod(bdflag(1)));
           aux=nodes(1,:);
           nodes(1,:)=nodes(m,:);
           nodes(m,:)=aux;
        end
        n=Old_Mesh.nAdjElements(i);
        for k=1:(n-1)
           l=k+1;
           while (length(setdiff(nodes(k,:),nodes(l,:))) ~= 1)
               l=l+1;
           end
           aux=nodes(l,:);
           nodes(l,:)=nodes(k+1,:);
           nodes(k+1,:)=aux;
           AdjNodes(i,k+1)=setdiff(intersect(nodes(k,:),nodes(k+1,:)),i);
        end
        if (~isempty(bdflag))
           AdjNodes(i,1)=u_nod(bdflag(1));
           AdjNodes(i,n+1)=u_nod(bdflag(2));
        else
           AdjNodes(i,1)=setdiff(intersect(nodes(1,:),nodes(n,:)),i);
        end
       end
   end
    
   
   
   % Assign output arguments
   
   New_Mesh = deal(Old_Mesh);
   New_Mesh.AdjNodes = AdjNodes;
   
return
