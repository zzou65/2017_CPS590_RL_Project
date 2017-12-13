function Mesh = smooth(Mesh,FixedPos)
% SMOOTH Laplacian smoothing.
%
%   MESH = SMOOTH(MESH,FIXEDPOS) smoothes the struct MESH using the Laplacian
%   smoothing algorithm.
%
%   The array FIXEDPOS specifies the fixed vertices of the mesh:
%    0 Vertex position will be moved during smoothing process.
%    1 Vertex position does not move during smoothing process.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements
%                of the mesh. 
%
%   Example:
%
%   Mesh = smooth(Mesh,FixedPos);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  nCoordinates = size(Mesh.Coordinates,1);
  [nElements,nVert] = size(Mesh.Elements); 
  
  % Initialize constants
  
  MAX_NODES = 10;           % Maximum number of adjacent nodes per patch
  MAX_EDGES = 3*nElements;  % Maximum number of edges
  TOL = 0.001;              % Mesh movement tolerance
  MAX_SWEEPS = 100;         % Maximum number of mesh sweeps
  
  % Build list of adjacent vertices and elements
  
  AdjNodes = zeros(nCoordinates,MAX_NODES);
  nAdjNodes = zeros(nCoordinates,1);
  if(nVert == 3)
      
    % Triangular elements  
      
    for i = 1:nElements
      for j = 1:3
        
        % Check wheter vertex 1 is a fixed position 
        
        idx_1 = Mesh.Elements(i,j);
        idx_2 = Mesh.Elements(i,rem(j+3,3)+1);
        idx_3 = Mesh.Elements(i,rem(j+4,3)+1);
      
        % Check for vertices 2 and 3 in the set of adjacent nodes of vertex 1
      
        sz = nAdjNodes(idx_1);
        tmp = AdjNodes(idx_1,:);
        occ_1 = 0;
        occ_2 = 0;
        k = 1;
        while((~occ_1 || ~occ_2) && k <= sz)
          occ_1 = max(occ_1,idx_2 == tmp(k));
          occ_2 = max(occ_2,idx_3 == tmp(k));
          k = k+1;
        end
        if(~occ_1)
          sz = sz+1;
          tmp(sz) = idx_2;    
        end
        if(~occ_2)
          sz = sz+1; 
          tmp(sz) = idx_3;
        end
   
        nAdjNodes(idx_1) = sz;
        AdjNodes(idx_1,:) = tmp;   
      end
    end
  else
  
    % Quadrilateral elements  
      
     for i = 1:nElements
      for j = 1:4
        
        % Check wheter vertex 1 is a fixed position
        
        idx_1 = Mesh.Elements(i,j);
        idx_2 = Mesh.Elements(i,rem(j+4,4)+1);
        idx_4 = Mesh.Elements(i,rem(j+6,4)+1);
      
        % Check for vertices 2 and 4 in the set of adjacent nodes of vertex 1
      
        sz = nAdjNodes(idx_1);
        tmp = AdjNodes(idx_1,:);
        occ_1 = 0;
        occ_2 = 0;
        k = 1;
        while((~occ_1 || ~occ_2) && k <= sz)
          occ_1 = max(occ_1,idx_2 == tmp(k));
          occ_2 = max(occ_2,idx_4 == tmp(k));
          k = k+1;
        end
        if(~occ_1)
          sz = sz+1;
          tmp(sz) = idx_2;    
        end
        if(~occ_2)
          sz = sz+1; 
          tmp(sz) = idx_4;
        end
           
        nAdjNodes(idx_1) = sz;
        AdjNodes(idx_1,:) = tmp;      
      end
    end    
  end
  
  % Squeeze adjacency lists
  
  MAX_NODES = max(nAdjNodes);
  AdjNodes = AdjNodes(:,1:MAX_NODES);
  
  % Start laplacian smoothing algorithm
  
  delta = inf*ones(nCoordinates,1);
  iter = 0;
  while(iter < MAX_SWEEPS && max(delta) > TOL)
    for i = 1:nCoordinates
      if(~FixedPos(i))  
        x_old = Mesh.Coordinates(i,:);      
        Mesh.Coordinates(i,:) = 1/nAdjNodes(i)*sum(Mesh.Coordinates(AdjNodes(i,1:nAdjNodes(i)),:),1);                     
        delta(i) =  sqrt(sum((x_old-Mesh.Coordinates(i,:)).^2));
      else      
        delta(i) = 0;
      end 
    end
    iter = iter+1;
  end      
  
return