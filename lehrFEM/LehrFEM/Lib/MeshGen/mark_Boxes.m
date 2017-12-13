function Mesh = mark_Boxes(Mesh,varargin)
% MARK_BOXES Mark elements of a mesh.
%
%   MESH = MARK_BOXES(MESH,BOX,FLAG) marks the elements of the struct MESH 
%   according to the location of their barycenters with the integer FLAG.
%   If the barycenter of an element is contained inside the box BOX its
%   element flag will be set to FLAG.
%
%   A box is defined by a 2-by-2 matrix whose first row represents the
%   starting point and the second box represents the end point.
%
%   MESH = MARK_BOXES(MESH,BOX1,FLAG1,BOX2,FLAG2,BOX3,FLAG3)
%
%   Example:
%
%   Mesh = mark_Boxes(Mesh,[0 0; 1 1],1);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nBoxes = (nargin()-1)/2;
  
  % Add element flags to the mesh
  
  Mesh.ElemFlag = zeros(nElements,1);
  
  % Compute barycenters of the mesh
  
  Xbar = (Mesh.Coordinates(Mesh.Elements(:,1),:) + ...
          Mesh.Coordinates(Mesh.Elements(:,2),:) + ...
          Mesh.Coordinates(Mesh.Elements(:,3),:))/3;
      
  % Extract boxes and flags
  
  Boxes = cell(nBoxes,1);
  Flags = zeros(nBoxes,1);
  j = 1;
  while(j <= nBoxes)
    Boxes{j} = varargin{2*j-1};
    Flags(j) = varargin{2*j};
    j = j+1;
  end
  
  % Assign element flags    
      
  for i = 1:nElements
    
    % Check wheter barycenter is contained inside current box    

    for j = 1:nBoxes
      if(isContained(Xbar(i,:),Boxes{j}(1,:),Boxes{j}(2,:)))
        Mesh.ElemFlag(i) = Flags(j);  
      end
    end
      
  end
  
return

% Subroutine: Checks wheter or not vertex is contained inside a box

function chk = isContained(vert,x0,x1)
% ISCONTAINED Checks wheter or not a vertex is contained inside a box.
%
%   CHK = ISCONTAINED(VERT,X0,X1) checks wheter the vertex VERT is
%   contained isnside the box with starting point X0 and final point X1.
%
%   Example:
%
%   isContained([0, 1],[0 0],[1 1]);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  chk = 0;
  if(vert(1) > x0(1) && vert(1) < x1(1))
    if(vert(2) > x0(2) && vert(2) < x1(2))
      chk = 1;   
    end
  end    
      
return