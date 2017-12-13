function Mesh = genMesh_PML(L,PML,varargin)
% GENMESH_PML Generate 2D PML Mesh.
%
%   MESH = GENMESH_PML(L,PML) generates a 2D PML mesh with on the
%   computational domain [-L,L] x [-L,L] with a PML of length PML.
%
%   MESH = GENMESH_PML(L,PML,HMAX) generates a 2D PML mesh with maximal
%   edge length less or equal to HMAX.
%
%   Example:
%
%   Mesh = genMesh(0.75,1,0.2);
%
%   See also orient_Elems.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute maximal edge length

  if(nargin < 3)
    hmax = 2*(L+PML)/10;
  else
    hmax = varargin{1};  
  end
  
  % Generate vertex coordinates

  h = PML/2;
  while(h > hmax)
    h = h/2;
  end
  xl = -(PML+L):h:-L;
    
  h = L;
  while(h > hmax)
    h = h/2;
  end
  xm = -L:h:L;
  
  h = PML/2;
  while(h > hmax)
    h = h/2;  
  end
  xr = L:h:(L+PML);
  
  xx = [xl ...
        xm(2:end) ...
        xr(2:end)];
  [x,y] = meshgrid(xx,xx);
  
  Mesh.Coordinates = [x(:) y(:)];
  
  % Generate Delaunay triangulation

  Mesh.Elements = delaunayn(Mesh.Coordinates);
  
  % Adjust element orientations
  
  Mesh = orient_Elems(Mesh);
  
return