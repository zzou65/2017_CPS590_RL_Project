function Mesh = TProd_Mesh(x,y)
%TPROD_MESH Create Tensor-Product mesh
%
%   MESH = TPROD_MESH(X,Y) initializes a tensor-product mesh with vertices
%   at (X(i),Y(j)).  The vectors X and Y must be sorted in ascending order.
%
%   MESH = TPROD_MESH(X) is short for MESH = TPROD_MESH(X,X);
%
%   Example :
%
%     Mesh = TProd_Mesh(0:1/n:1);

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Handle input arguments
  if(nargin < 2)
    y = x;
  end
  
  % Create coordinates
  [X,Y] = meshgrid(x,y);
  Mesh.Coordinates = [X(:),Y(:)];
  
  % Create elements
  nX = length(x);
  nY = length(y);
  nElem = (nX-1)*(nY-1);
  Mesh.Elements = zeros(nElem,4);
  for i=1:nX-1
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),1) = (i-1)*nY + (1:nY-1);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),4) = (i-1)*nY + (2:nY);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),2) = i*nY + (1:nY-1);
    Mesh.Elements((i-1)*(nY-1)+(1:nY-1),3) = i*nY + (2:nY);
  end

return