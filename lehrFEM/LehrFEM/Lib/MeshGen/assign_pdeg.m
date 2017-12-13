function [EDofs,CDofs,ElemDeg] = assign_pdeg(Mesh,CNodes,pmax)
% ASSIGN_PDEG Assign polynomial degrees.
%
%   [EDOFS,CDOFS,ELEMDEG] = ASSIGN_PDEG(MESH,CNODES,PMAX) assigns the polynomial
%   degrees to the edges and elements of an adaptively refined mesh with
%   the corner nodes CNODES. PMAX specifies the maximum polynomial degree
%   on the mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.
%
%   Example:
%
%   [EDofs,CDofs,ElemDeg] = assign_pdeg(Mesh,CNodes,pmax);
%
%   See also add_Patches.

%   Copyright 2006-2006 Patrick Meury & Mengyu Wang & Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCNodes = size(CNodes,1);           % Number of corner nodes
  nEdges = size(Mesh.Edges,1);        % Number of edges
  nElements = size(Mesh.Elements,1);  % Number of elements
   
  % Add patches to current mesh
  
  Mesh = add_Patches(Mesh);
  
  % Assign degree of polynomials to elements
  
  ElemDeg = -ones(nElements,nCNodes); 
  
  % compute the polynomial degree for every element with respect to every
  % corner node
  for i = 1:nCNodes 
    
    p = 1;
    NodeList = CNodes(i);
    
    % Start Loop over all nodes of the mesh
    
    while(1)
      
      NewNodes = [ ];
      
      % iterate over all nodes in NodeList
      while(~isempty(NodeList)) 
      
        % Assign polynomial degree to adjacent elements  
        
        node = NodeList(1);
        for j = 1:Mesh.nAdjElements(node)
          elem = Mesh.AdjElements(node,j);
          
          % if the element has not yet been processed assign the polynomial
          % degree p
          if(ElemDeg(elem,i) == -1)
            ElemDeg(elem,i) = p;  
            
            % add the new nodes of the current element to NewNodes
            NewNodes = [NewNodes setdiff(Mesh.Elements(elem,:),node)];
          end
        end
      
        % Remove current node from list 
        
        NodeList = NodeList(2:end);
       
      end
      
      % Increase polynomial degree
      
      p = min(p+1,pmax);
      
      % Break Loop if all nodes have been checked
      
      NodeList = NewNodes;
      if(isempty(NodeList))
        break;  
      end
      
    end
  end
  
  ElemDeg = min(ElemDeg,[],2);
  
  % assign inner degrees of freedom
  
  CDofs = zeros(nElements,1);
  loc = ElemDeg>2;
  CDofsGl = (ElemDeg-1).*(ElemDeg-2)/2;
  CDofs(loc) = CDofsGl(loc);
  
 
  % Assign degree of polynomials to edges
  % computation according to maximum of the two adjacent elements
  
  EDofs = zeros(nEdges,1);
  for i = 1:nEdges      
    if(Mesh.Edge2Elem(i,1))
      p_left = ElemDeg(Mesh.Edge2Elem(i,1))-1; 
    else
      p_left = 0;
    end  
    if(Mesh.Edge2Elem(i,2))
      p_right = ElemDeg(Mesh.Edge2Elem(i,2))-1;  
    else
      p_right = 0;
    end
    EDofs(i) = max(p_left,p_right);
  end
  
  
return  